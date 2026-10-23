"""Read mutation residue templates from the supplied simulation force field."""

import copy
from contextlib import contextmanager
import io
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import warnings
from types import SimpleNamespace

import parmed


def mutation_atom_names(residue, mutant):
    """Canonical names for matching AMBER/GROMACS/CHARMM residue templates."""
    names = {atom.name for atom in residue.atoms}
    aliases = {'HN': 'H', 'HT1': 'H1', 'HT2': 'H2', 'HT3': 'H3',
               'OT1': 'O', 'OT2': 'OXT', 'O1': 'O', 'O2': 'OXT',
               'OC1': 'O', 'OC2': 'OXT'}
    if mutant == 'GLY' and 'HA1' not in names:
        aliases['HA3'] = 'HA1'
    return [aliases.get(atom.name, atom.name) for atom in residue.atoms]


class MutationForceField:
    """Residue charges/types plus the actual nonbonded parameter definitions."""

    def __init__(self, atom_types=None, combining_rule='lorentz', pair_types=None):
        self.atom_types = copy.deepcopy(atom_types or {})
        self.pair_types = dict(pair_types or {})
        for name, atom_type in self.atom_types.items():
            pair = self.pair_types.get((name, name))
            if pair is not None:
                atom_type.rmin_14 = pair.rmin / 2
                atom_type.epsilon_14 = pair.epsilon
        self.combining_rule = combining_rule
        self.records = []
        self.atom_masses = {}
        self.sources = []

    @staticmethod
    def gromacs_library_dirs(command=None):
        """Find data for the selected GROMACS, including relocated installations."""
        roots = []
        if os.environ.get('GMXLIB'):
            roots.extend(Path(p).expanduser() for p in os.environ['GMXLIB'].split(os.pathsep) if p)
        if command:
            executable = str(command[0] if isinstance(command, (list, tuple)) else command)
            executable = Path(shutil.which(executable) or executable).expanduser().absolute()
            # The reported prefix handles wrappers whose own location is unrelated
            # to the installation. Inspect both output streams across versions.
            try:
                result = subprocess.run([str(executable), '--version'], capture_output=True,
                                        text=True, check=False, timeout=5)
                output = result.stdout + '\n' + result.stderr
                match = re.search(r'^\s*Data prefix:\s*(.+?)\s*$', output, re.MULTILINE)
                if match:
                    prefix = Path(re.sub(r'\s+\(source tree\)$', '', match[1])).expanduser()
                    roots.extend((prefix / 'share/gromacs/top', prefix / 'share/top'))
            except (OSError, subprocess.TimeoutExpired):
                pass
            # Also support symlinks, older executables and unavailable version output.
            for binary in (executable.resolve(), executable):
                prefix = binary.parent.parent
                roots.extend((prefix / 'share/gromacs/top', prefix / 'share/top'))
        roots.append(Path(sys.prefix) / 'share/gromacs/top')
        topdir = getattr(parmed.gromacs, 'GROMACS_TOPDIR', None)
        if topdir:
            roots.append(Path(topdir))
        return tuple(dict.fromkeys(path.resolve() for path in roots if path.is_dir()))

    @staticmethod
    @contextmanager
    def gromacs_library(directories):
        """Use the selected installation for ParmEd's global include directory."""
        if not directories:
            yield
            return
        previous = parmed.gromacs.GROMACS_TOPDIR
        try:
            parmed.gromacs.GROMACS_TOPDIR = str(directories[0])
            yield
        finally:
            parmed.gromacs.GROMACS_TOPDIR = previous

    @classmethod
    def candidate_sources(cls, topology, gromacs_dirs=()):
        types = {a.type for a in topology.atoms}
        gmxroots = gromacs_dirs or cls.gromacs_library_dirs()
        if getattr(topology, 'chamber', False) or {'NH1', 'CT1'} <= types:
            family = 'charmm'
        elif any(name.startswith('opls') for name in types):
            family = 'opls'
        else:
            family = 'amber'
        sources = []
        for root in gmxroots:
            for path in sorted(Path(root).glob('*.ff')):
                known_family = next((f for f in ('amber', 'charmm', 'opls', 'gromos')
                                     if path.name.lower().startswith(f)), None)
                if known_family in (None, family) and (path / 'forcefield.itp').is_file():
                    sources.append(path)
        if family != 'amber':
            return sources
        for name in ('leaprc.protein.ff14SB', 'leaprc.protein.ff03.r1', 'leaprc.ff99SB',
                     'leaprc.protein.ff19SB', 'leaprc.protein.ff15ipq'):
            try:
                sources.append(MutationForceField._amber_file(name, Path.cwd()))
            except FileNotFoundError:
                pass
        return sources

    def fit(self, topology):
        """Rank compatible sources using residue, hydrogen, net-charge and LJ patterns."""
        residuals, hydrogen, net = [], [], []
        for residue in topology.residues:
            base = residue.name[1:] if len(residue.name) == 4 and residue.name[0] in 'NC' else residue.name
            actual_names = mutation_atom_names(residue, base)
            actual = dict(zip(actual_names, residue.atoms))
            choices = []
            for name, atoms, _ in self.records:
                if name not in (base, 'N' + base, 'C' + base):
                    continue
                dummy = SimpleNamespace(atoms=[SimpleNamespace(name=n) for n in atoms])
                keys = mutation_atom_names(dummy, base)
                pattern = dict(zip(keys, atoms.values()))
                if pattern.keys() == actual.keys():
                    dq = [actual[key].charge - charge for key, (_, charge) in pattern.items()]
                    dh = [actual[key].charge - charge for key, (_, charge) in pattern.items()
                          if actual[key].atomic_number == 1]
                    choices.append((sum(q*q for q in dq), dq, dh))
            if choices:
                _, dq, dh = min(choices, key=lambda item: item[0])
                residuals.extend(dq)
                hydrogen.extend(dh)
                net.append(sum(dq))
        lj = []
        for atom in topology.atoms:
            ref = self.atom_types.get(atom.type)
            if ref is not None and ref.rmin is not None and ref.epsilon is not None:
                lj.extend(((atom.rmin - ref.rmin) / max(ref.rmin, .1),
                           (atom.epsilon - ref.epsilon) / max(ref.epsilon, .01)))
        rms = lambda values: math.sqrt(sum(x*x for x in values) / len(values)) if values else 0.
        charge_rms = rms(residuals)
        score = charge_rms + rms(hydrogen) + rms(net) + .05*rms(lj)
        if not residuals:
            # Sparse/unusual residues provide weaker evidence; keep this visible.
            score += 1.
        return score, charge_rms, rms(hydrogen), len(residuals)

    @classmethod
    def infer(cls, topology, mutant, names, gromacs_parameters=None, gromacs_dirs=()):
        candidates = []
        for source in dict.fromkeys(cls.candidate_sources(topology, gromacs_dirs)):
            try:
                with warnings.catch_warnings(), cls.gromacs_library(gromacs_dirs):
                    warnings.simplefilter('ignore', parmed.exceptions.ParameterWarning)
                    warnings.simplefilter('ignore', parmed.exceptions.AmberWarning)
                    library = cls.load(str(source))
                matches = [r for r in library.topology(mutant).residues
                           if frozenset(mutation_atom_names(r, mutant)) == names]
                if matches:
                    candidates.append((library.fit(topology), library, source))
            except (OSError, ValueError, parmed.exceptions.ParameterError):
                continue
        if not candidates:
            return None, None
        fit, library, source = min(candidates, key=lambda item: item[0][0])
        # Retain known source-topology LJ parameters. The inferred library fills
        # missing residue charges and types, rather than replacing known LJ data.
        actual = dict(getattr(gromacs_parameters, 'atom_types', {}) or {})
        for atom in topology.atoms:
            actual.setdefault(atom.type, atom.atom_type)
        # Compiled component topologies may retain unused types (e.g. GLY HB2)
        # whose pair tables are still authoritative even with no such atom.
        for typename, index in getattr(topology, 'LJ_types', {}).items():
            if typename in actual:
                continue
            prototype = library.atom_types.get(typename)
            if prototype is None:
                element = 1 if typename.startswith('H') else 6 if typename.startswith('C') else None
                if element is None:
                    continue
                prototype = parmed.AtomType(typename, None, 1.008 if element == 1 else 12.011, element)
            prototype = copy.deepcopy(prototype)
            prototype.set_lj_params(topology.LJ_depth[index-1], topology.LJ_radius[index-1],
                                    topology.LJ_14_depth[index-1] if topology.chamber else None,
                                    topology.LJ_14_radius[index-1] if topology.chamber else None)
            actual[typename] = prototype
        for typename in library.atom_types.keys() & actual.keys():
            library.atom_types[typename] = copy.deepcopy(actual[typename])
        if 'HB1' in actual and ('NH1' in actual or getattr(topology, 'chamber', False)):
            # CHARMM36 distinguishes alpha CH, alpha CH2, and methyl hydrogens.
            # Older compatible charge libraries use HA/HB for multiple roles.
            # Resolve the modern roles from existing/source LJ definitions.
            for modern, older in (('HB1', 'HB'), ('HB2', 'HA2'), ('HA3', 'HA3')):
                prototype = actual.get(modern) or library.atom_types.get(older)
                if prototype is not None:
                    library.atom_types[modern] = copy.deepcopy(prototype)
                    for (a, b), pair in list(library.pair_types.items()):
                        if a == older or b == older:
                            library.pair_types[(modern if a == older else a,
                                                modern if b == older else b)] = pair
            records = []
            for name, atoms, origin in library.records:
                if name in ('ALA', 'NALA', 'CALA', 'GLY', 'NGLY', 'CGLY'):
                    atoms = dict(atoms)
                    for key, (typ, charge) in atoms.items():
                        desired = ('HB2' if name.endswith('GLY') and key.startswith('HA') else
                                   'HA3' if name.endswith('ALA') and key.startswith('HB') else
                                   'HB1' if name.endswith('ALA') and key == 'HA' else typ)
                        if desired in library.atom_types:
                            atoms[key] = (desired, charge)
                records.append((name, atoms, origin))
            library.records = records
        return library, dict(source=str(source), score=fit[0], charge_rms=fit[1],
                             hydrogen_rms=fit[2], matched_atoms=fit[3],
                             molecule_charge=sum(a.charge for a in topology.atoms))

    @staticmethod
    def _amber_file(name, parent):
        roots = [Path(parent)]
        amberhome = os.environ.get('AMBERHOME') or parmed.amber.AMBERHOME
        if amberhome:
            leap = Path(amberhome) / 'dat/leap'
            roots.extend(leap / directory for directory in ('cmd', 'cmd/oldff', 'lib', 'lib/oldff', 'parm'))
        for root in roots:
            candidate = root / name
            if candidate.is_file():
                return candidate.resolve()
        raise FileNotFoundError(f'Cannot find force-field file {name}')

    def _read_leaprc(self, path):
        from parmed.amber import AmberParameterSet

        def expand(filename, stack=()):
            filename = Path(filename).resolve()
            if filename in stack:
                raise ValueError(f'Cyclic force-field source reference: {filename}')
            lines = []
            for line in filename.read_text().splitlines():
                words = shlex.split(line.split('#', 1)[0])
                if words and words[0].lower() == 'source':
                    child = self._amber_file(words[1], filename.parent)
                    lines.extend(expand(child, (*stack, filename)))
                elif words and words[0].lower() in ('loadoff', 'loadamberparams'):
                    source = self._amber_file(words[1], filename.parent)
                    lines.append(f'{words[0]} "{source}"')
                else:
                    lines.append(line)
            self.sources.append(str(filename))
            return lines

        # One parameter set preserves LEaP's load order and later overrides.
        # Resolve relative paths without changing the process working directory.
        params = AmberParameterSet.from_leaprc(io.StringIO('\n'.join(expand(path))), search_oldff=True)
        self.atom_types.update(params.atom_types)
        for name, residue in params.residues.items():
            if not hasattr(residue, 'atoms'):
                continue  # Multi-residue library units are not residue templates.
            self.records.append((name, {a.name: (a.type, a.charge) for a in residue.atoms}, str(path)))

    def _read_rtp(self, path):
        residues = {}
        current = section = None
        subsections = {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                       'exclusions', 'cmap', 'bondedtypes'}
        for raw in path.read_text().splitlines():
            line = raw.split(';', 1)[0].strip()
            if not line or line.startswith('#'):
                continue
            header = re.fullmatch(r'\[\s*([^]]+?)\s*\]', line)
            if header:
                section = header[1].strip()
                if section.lower() not in subsections:
                    current = section
                    residues.setdefault(current, {})
                continue
            if current and section == 'atoms':
                words = line.split()
                residues[current][words[0]] = (words[1], float(words[2]))
        mass_table = {}
        atp = path.parent / 'atomtypes.atp'
        if atp.is_file():
            for raw in atp.read_text().splitlines():
                words = raw.split(';', 1)[0].split()
                if len(words) >= 2 and not words[0].startswith('#'):
                    mass_table[words[0]] = float(words[1])
        for name, atoms in residues.items():
            if atoms:
                self.records.append((name, atoms, str(path)))
                self.atom_masses[str(path), name] = {
                    atom: mass_table[typ] for atom, (typ, _) in atoms.items() if typ in mass_table}
        self.sources.append(str(path))
        for suffix, prefix in (('.n.tdb', 'N'), ('.c.tdb', 'C')):
            tdb = path.with_suffix(suffix)
            if tdb.is_file():
                self._read_tdb(tdb, residues, prefix, str(path))

    def _read_tdb(self, path, residues, prefix, rtp_source):
        patches = {}
        current = section = None
        for raw in path.read_text().splitlines():
            line = raw.split(';', 1)[0].strip()
            if not line or line.startswith('#'):
                continue
            header = re.fullmatch(r'\[\s*([^]]+?)\s*\]', line)
            if header:
                section = header[1].strip().lower()
                if section not in ('replace', 'add', 'delete', 'bonds', 'angles', 'dihedrals', 'impropers'):
                    current = header[1].strip()
                    patches[current] = {}
                continue
            if current:
                patches[current].setdefault(section, []).append(line.split())
        for target in ('ALA', 'GLY'):
            if target not in residues:
                continue
            for patchname, patch in patches.items():
                if '-' in patchname and patchname.split('-', 1)[0] in ('ALA', 'GLY', 'PRO'):
                    if not patchname.startswith(target + '-'):
                        continue
                elif target + '-' + patchname in patches:
                    continue
                atoms = dict(residues[target])
                masses = dict(self.atom_masses.get((rtp_source, target), {}))
                try:
                    for row in patch.get('delete', []):
                        del atoms[row[0]]
                    for row in patch.get('replace', []):
                        old = row[0]
                        if old not in atoms and not (prefix == 'C' and old == 'OXT'):
                            raise KeyError(old)
                        atoms.pop(old, None)
                        masses.pop(old, None)
                        name = row[1] if len(row) == 5 else old
                        atoms[name] = (row[-3], float(row[-1]))
                        masses[name] = float(row[-2])
                    adds = iter(patch.get('add', []))
                    for row in adds:
                        count, _, name = row[:3]
                        parameters = next(adds)
                        for i in range(int(count)):
                            atoms[name + (str(i + 1) if int(count) > 1 else '')] = (
                                parameters[0], float(parameters[2]))
                            masses[name + (str(i + 1) if int(count) > 1 else '')] = float(parameters[1])
                except (KeyError, ValueError, StopIteration):
                    # A patch for another residue/terminal chemistry is not a match.
                    continue
                if atoms != residues[target] and abs(sum(q for _, q in atoms.values()) -
                                                     round(sum(q for _, q in atoms.values()))) < 1e-4:
                    source = f'{path}:{patchname}'
                    self.records.append((prefix + target, atoms, source))
                    self.atom_masses[source, prefix + target] = masses

    def _read_charmm(self, paths, use_parameters=True):
        from parmed.charmm import CharmmParameterSet

        params = CharmmParameterSet()
        for path in paths:
            if path.suffix.lower() == '.rtf':
                with warnings.catch_warnings():
                    # Donor/acceptor annotations do not enter charges or LJ.
                    warnings.filterwarnings('ignore', r'WARNING: Ignoring "DELETE (ACCE|DONO) ',
                                            parmed.exceptions.ParameterWarning)
                    params.read_topology_file(str(path))
        if use_parameters:
            for path in paths:
                if path.suffix.lower() == '.prm':
                    params.read_parameter_file(str(path))
        for path in paths:
            if path.suffix.lower() == '.str':
                params.read_stream_file(str(path))
        # The GROMACS include tree, when present, defines the effective LJ data.
        self.atom_types = {**params.atom_types, **self.atom_types}
        for residue in (*params.residues.values(), *params.patches.values()):
            for atom in residue.atoms:
                if atom.type in self.atom_types:
                    atom.atomic_number = self.atom_types[atom.type].atomic_number
                    atom.mass = self.atom_types[atom.type].mass
        source = ', '.join(map(str, paths))
        for name, residue in params.residues.items():
            self.records.append((name, {a.name: (a.type, a.charge) for a in residue.atoms}, source))
            if name not in ('ALA', 'GLY'):
                continue
            for prefix, patchname in (('N', 'GLYP' if name == 'GLY' else 'NTER'), ('C', 'CTER')):
                if patchname not in params.patches:
                    continue
                patched = residue.apply_patch(params.patches[patchname])
                self.records.append((prefix + name,
                                     {a.name: (a.type, a.charge) for a in patched.atoms}, source))
        self.sources.extend(map(str, paths))

    @classmethod
    def load(cls, explicit='', topology_file=None, includes=(), gromacs_parameters=None,
             combining_rule='lorentz'):
        library = cls(getattr(gromacs_parameters, 'atom_types', {}), combining_rule,
                      getattr(gromacs_parameters, 'pair_types', {}))
        if explicit:
            path = Path(explicit).expanduser()
            if not path.exists():
                path = cls._amber_file(explicit, Path.cwd())
            roots = {path} if path.is_dir() else {path.parent}
            files = [path] if path.is_file() else []
            single_file = path.is_file()
        else:
            single_file = False
            files = []
            roots = {Path(p).resolve().parent for p in includes}
            if topology_file:
                parent = Path(topology_file).resolve().parent
                roots.update((parent, parent / 'toppar'))
                if parent.name.lower() in ('gromacs', 'amber'):
                    roots.add(parent.parent / 'toppar')
            if not roots:
                return library
        for directory in sorted(roots):
            if not directory.is_dir():
                continue
            if not single_file:
                files.extend(sorted(directory.glob('*.rtp')))
                files.extend(sorted(directory.glob('top*prot*.rtf')))
                if gromacs_parameters is None:
                    files.extend(sorted(directory.glob('par*prot*.prm')))
        files = list(dict.fromkeys(files))
        # An explicitly selected GROMACS force-field directory supplies its own
        # atom types even for a native AMBER input.
        if explicit and Path(explicit).is_dir() and (Path(explicit) / 'forcefield.itp').is_file():
            gmx = parmed.gromacs.GromacsTopologyFile(str(Path(explicit) / 'forcefield.itp'), parametrize=False)
            parsed = cls(gmx.parameterset.atom_types, gmx.combining_rule, gmx.parameterset.pair_types)
            library.atom_types = parsed.atom_types
            library.pair_types = parsed.pair_types
            library.combining_rule = gmx.combining_rule
        charmm = [p for p in files if p.suffix.lower() in ('.rtf', '.prm', '.str')]
        if charmm:
            library._read_charmm(charmm, use_parameters=gromacs_parameters is None)
        for path in files:
            if path.suffix.lower() == '.rtp':
                library._read_rtp(path)
            elif path.name.startswith('leaprc'):
                library._read_leaprc(path)
        return library

    def validate_topology(self, topology):
        """Reject a source that contradicts identifiable residues/LJ types."""
        for residue in topology.residues:
            actual = {a.name: a.charge for a in residue.atoms}
            candidates = [atoms for name, atoms, _ in self.records
                          if name in (residue.name, 'N' + residue.name, 'C' + residue.name)
                          and atoms.keys() == actual.keys()]
            if candidates and not any(all(abs(actual[name] - q) <= 1e-4
                                           for name, (_, q) in pattern.items()) for pattern in candidates):
                raise ValueError(f'Mutation force field does not match the source topology charges '
                                 f'for residue {residue.name} {residue.idx + 1}')
        for atom in topology.atoms:
            reference = self.atom_types.get(atom.type)
            if reference is not None and reference.rmin is not None and reference.epsilon is not None:
                if not (math.isclose(atom.rmin, reference.rmin, rel_tol=1e-5, abs_tol=1e-6) and
                        math.isclose(atom.epsilon, reference.epsilon, rel_tol=1e-5, abs_tol=1e-6)):
                    raise ValueError(f'Mutation force field does not match the source topology LJ '
                                     f'parameters for atom type {atom.type}')

    def topology(self, mutant):
        residues = []
        for name, atoms, source in self.records:
            if name not in (mutant, 'N' + mutant, 'C' + mutant):
                continue
            result = []
            for atomname, (typename, charge) in atoms.items():
                atom_type = self.atom_types.get(typename)
                if atom_type is None or atom_type.rmin is None or atom_type.epsilon is None:
                    break
                atom = parmed.Atom(name=atomname, type=typename, charge=charge,
                                   mass=self.atom_masses.get((source, name), {}).get(atomname, atom_type.mass),
                                   atomic_number=atom_type.atomic_number)
                atom.atom_type = atom_type
                result.append(atom)
            else:
                residues.append(SimpleNamespace(name=name, atoms=result, source=source))
        return SimpleNamespace(residues=residues)


def missing_lj14_pairs(parm, original, residue):
    """Newly used Chamber 1-4 type pairs absent from a sparse source table."""
    if not parm.chamber:
        return []
    active = {tuple(sorted((d.atom1.nb_idx, d.atom4.nb_idx))) for d in original.dihedrals
              if not d.ignore_end and not d.improper}
    missing = {}
    for dihedral in parm.dihedrals:
        if dihedral.ignore_end or dihedral.improper:
            continue
        first, second = dihedral.atom1, dihedral.atom4
        if first.residue is not residue and second.residue is not residue:
            continue
        pair = tuple(sorted((first.nb_idx, second.nb_idx)))
        if pair in active:
            continue
        slot = parm.parm_data['NONBONDED_PARM_INDEX'][
            parm.ptr('ntypes') * (first.nb_idx - 1) + second.nb_idx - 1] - 1
        if slot >= 0 and all(parm.parm_data[flag][slot] == 0 for flag in
                             ('LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF')):
            missing[slot] = (first, second)
    return [(slot, *atoms) for slot, atoms in missing.items()]


def fill_lj14_pairs(parm, pairs, library):
    """Fill only newly needed pairs, retaining all previously active interactions."""
    estimated = False
    for slot, first, second in pairs:
        explicit = (library.pair_types.get((first.type, second.type)) or
                    library.pair_types.get((second.type, first.type)))
        types = [library.atom_types.get(atom.type) for atom in (first, second)]
        nbfix = types[0].nbfix.get(second.type) if types[0] is not None else None
        if nbfix is None and types[1] is not None:
            nbfix = types[1].nbfix.get(first.type)
        if explicit is not None:
            radius, epsilon = explicit.rmin, explicit.epsilon
        elif nbfix:
            radius, epsilon = nbfix[2:4]
        else:
            radii, depths = [], []
            for atom, atom_type in zip((first, second), types):
                self_pair = library.pair_types.get((atom.type, atom.type))
                if self_pair is not None:
                    r, e = self_pair.rmin / 2, self_pair.epsilon
                else:
                    if atom_type is None:
                        atom_type = atom.atom_type
                        estimated = True
                    r, e = atom_type.rmin_14, atom_type.epsilon_14
                    if r is None:
                        r = atom_type.rmin
                    if e is None:
                        e = atom_type.epsilon
                radii.append(r)
                depths.append(e)
            radius = (2 * math.sqrt(radii[0] * radii[1]) if library.combining_rule == 'geometric'
                      else sum(radii))
            epsilon = math.sqrt(depths[0] * depths[1])
        parm.parm_data['LENNARD_JONES_14_ACOEF'][slot] = epsilon * radius**12
        parm.parm_data['LENNARD_JONES_14_BCOEF'][slot] = 2 * epsilon * radius**6
    return estimated


def ensure_lj_type(parm, parameters, sample_atom, indices, atom_types=None,
                   combining_rule='lorentz', pair_types=None):
    """Resolve or add a type, preserving existing pairs and source NBFIX/1-4 data."""
    from parmed.tools.addljtype import AddLJType

    typename = parameters['type']
    available = indices.get(typename, set())
    if len(available) == 1:
        return next(iter(available))
    if available:
        raise ValueError(f'Ambiguous Lennard-Jones type {typename}: {len(available)} nonbonded indices')
    if any(parm.parm_data.get('LENNARD_JONES_CCOEF', [])):
        raise ValueError('Adding a missing mutation LJ type to a 12-6-4 topology requires explicit '
                         'source C-coefficients; they cannot be inferred from 12-6 parameters')
    atom_type = parameters['atom_type']
    atom_types = atom_types or {}
    pair_types = pair_types or {}
    flags = ['LENNARD_JONES_ACOEF', 'LENNARD_JONES_BCOEF']
    if parm.chamber:
        flags += ['LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF']

    def coefficients(othername, index):
        if othername in parameters.get('pairs', {}):
            return parameters['pairs'][othername]
        other = atom_types.get(othername)
        nbfix = atom_type.nbfix.get(othername)
        if nbfix is None and other is not None:
            nbfix = other.nbfix.get(typename)
        values = []
        for offset in range(2 if parm.chamber else 1):
            explicit14 = pair_types.get((typename, othername)) or pair_types.get((othername, typename))
            if offset and explicit14 is not None:
                radius, epsilon = explicit14.rmin, explicit14.epsilon
            elif nbfix:
                radius, epsilon = nbfix[2*offset:2*offset+2]
            else:
                radii = parm.LJ_14_radius if offset else parm.LJ_radius
                depths = parm.LJ_14_depth if offset else parm.LJ_depth
                own_r = atom_type.rmin_14 if offset else atom_type.rmin
                own_e = atom_type.epsilon_14 if offset else atom_type.epsilon
                own_r = atom_type.rmin if own_r is None else own_r
                own_e = atom_type.epsilon if own_e is None else own_e
                radius = (2 * math.sqrt(own_r * radii[index-1]) if combining_rule == 'geometric'
                          else own_r + radii[index-1])
                epsilon = math.sqrt(own_e * depths[index-1])
            values.extend((epsilon * radius**12, 2 * epsilon * radius**6))
        return tuple(values)

    def equal(a, b):
        return len(a) == len(b) and all(math.isclose(x, y, rel_tol=1e-7, abs_tol=1e-8)
                                       for x, y in zip(a, b))

    ntypes = parm.ptr('ntypes')
    rows = []
    for index in range(1, ntypes+1):
        aliases = [name for name, choices in indices.items() if index in choices] or ['']
        choices = [coefficients(name, index) for name in aliases]
        if not all(equal(choices[0], choice) for choice in choices[1:]):
            raise ValueError(f'Cannot add LJ type {typename}: existing index {index} combines atom types '
                             'with distinct source pair parameters')
        rows.append(choices[0])
    self_pair = list(parameters.get('pairs', {}).get(typename, ()))
    if not self_pair:
        nbfix = atom_type.nbfix.get(typename)
        for offset in range(2 if parm.chamber else 1):
            if offset and (typename, typename) in pair_types:
                radius = pair_types[typename, typename].rmin
                epsilon = pair_types[typename, typename].epsilon
            elif nbfix:
                radius, epsilon = nbfix[2*offset:2*offset+2]
            else:
                radius = atom_type.rmin_14 if offset else atom_type.rmin
                epsilon = atom_type.epsilon_14 if offset else atom_type.epsilon
                radius = 2 * (atom_type.rmin if radius is None else radius)
                epsilon = atom_type.epsilon if epsilon is None else epsilon
            self_pair.extend((epsilon * radius**12, 2 * epsilon * radius**6))

    # Two force-field names can legitimately share the same complete LJ row.
    # Match all interactions, including 1-4/NBFIX, before reusing an index.
    for candidate in range(1, ntypes+1):
        observed = []
        for other in range(1, ntypes+1):
            pair = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*(candidate-1)+other-1] - 1
            if pair < 0:
                break
            observed.append(tuple(parm.parm_data[flag][pair] for flag in flags))
        if len(observed) == ntypes and all(equal(a, b) for a, b in zip(observed, rows)) and \
                equal(observed[candidate-1], self_pair):
            indices[typename] = {candidate}
            parm.LJ_types[typename] = candidate
            return candidate

    # Plan/validate all rows first. AddLJType then extends the existing arrays;
    # only its newly added entries are replaced with source pair coefficients.
    parm.remake_parm()
    selection = [int(atom is sample_atom) for atom in parm.atoms]
    r14 = atom_type.rmin if atom_type.rmin_14 is None else atom_type.rmin_14
    e14 = atom_type.epsilon if atom_type.epsilon_14 is None else atom_type.epsilon_14
    AddLJType(parm, selection, atom_type.rmin, atom_type.epsilon, r14, e14)
    new_index = parm.ptr('ntypes')
    for other, values in enumerate([*rows, self_pair], start=1):
        pair = parm.parm_data['NONBONDED_PARM_INDEX'][new_index*(new_index-1)+other-1] - 1
        for flag, value in zip(flags, values):
            parm.parm_data[flag][pair] = value
    sample_atom.nb_idx = new_index
    parm.LJ_types[typename] = new_index
    indices[typename] = {new_index}
    return new_index
