from pathlib import Path
import re
import tempfile


INCLUDE_RE = re.compile(r'^(\s*#include\s+)([<"])([^>"]+)([>"])(.*)$')


def comment_gromacs_cmap(line, in_cmap=False):
    stripped = line.lstrip()
    active = bool(stripped) and not stripped.startswith(';')
    section = re.match(r'\[\s*([^\]]+?)\s*\]', stripped) if active else None

    if section:
        in_cmap = section.group(1).strip().lower() == 'cmap'

    if in_cmap:
        return line if stripped.startswith(';') else f';{line}', True, True

    return line, False, False


class GromacsTopologyPreprocessor:
    def __init__(self):
        self.created_files = []
        self._processed = {}
        self.cmap_found = False

    def preprocess(self, source, remove_solvent=False, solvent_ions=None):
        source = Path(source)
        if source in self._processed:
            return self._processed[source]

        with tempfile.NamedTemporaryFile(
            dir=source.parent, prefix=f'_temp_{source.stem}_', suffix=source.suffix,
            mode='w', delete=False
        ) as temp_file:
            temp_path = Path(temp_file.name)
            self.created_files.append(temp_path)
            self._processed[source] = temp_path

            molsect = False
            in_cmap = False
            with open(source) as input_file:
                for line in input_file:
                    line, in_cmap, cmap_line = comment_gromacs_cmap(line, in_cmap)
                    self.cmap_found = self.cmap_found or cmap_line
                    line = self._rewrite_local_include(line, source.parent)

                    if '[ molecules ]' in line:
                        molsect = True
                    if remove_solvent and molsect:
                        if not line.split():
                            continue
                        if line.split()[0].strip() in solvent_ions:
                            continue

                    temp_file.write(line)

        return temp_path

    def _rewrite_local_include(self, line, parent):
        line_end = '\n' if line.endswith('\n') else ''
        content = line[:-1] if line_end else line
        match = INCLUDE_RE.match(content)
        if not match:
            return line

        prefix, opener, include_name, closer, suffix = match.groups()
        include_path = Path(include_name)
        if include_path.is_absolute():
            return line

        resolved = parent / include_path
        if not resolved.exists():
            return line

        processed = self.preprocess(resolved)
        try:
            include_target = processed.relative_to(parent)
        except ValueError:
            include_target = processed
        return f'{prefix}{opener}{include_target.as_posix()}{closer}{suffix}{line_end}'
