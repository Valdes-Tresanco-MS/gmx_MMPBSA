import unittest
from pathlib import Path
import copy
import re
import tempfile

from GMXMMPBSA.input_parser import input_file


ROOT = Path(__file__).resolve().parents[1]


class LoggingDocumentationTest(unittest.TestCase):
    def test_logging_guide_documents_supported_modes_and_cluster_monitoring(self):
        guide = (ROOT / 'docs/logging.md').read_text()
        for text in ('gmx_MMPBSA.log', 'tail -f gmx_MMPBSA.log', '`auto`', '`rich`', '`classic`', '`plain`', '`none`'):
            with self.subTest(text=text):
                self.assertIn(text, guide)

    def test_logging_guide_is_in_mkdocs_navigation(self):
        nav = (ROOT / 'mkdocs.yml').read_text()
        self.assertIn('Logging and progress: logging.md', nav)

    def test_changelog_mentions_record_based_logging_changes(self):
        changelog = (ROOT / 'docs/changelog.md').read_text()
        self.assertIn('record-based warning/error totals', changelog)
        self.assertIn('not a stable machine-readable interface', changelog)


class CompatibilityMigrationDocumentationTest(unittest.TestCase):
    def test_170_migration_guide_covers_behavior_changes(self):
        guide = (ROOT / 'docs/compatibility.md').read_text()
        required = (
            '1.7.0 release',
            '1.6.5 environment',
            'Python `>=3.11,<3.13`',
            '`amber_MMPBSA`',
            'Composite alanine/glycine mutations',
            '`qh_entropy=1` is rejected',
            'full selected ensemble',
            '`Block SD`',
            '`Block SEM`',
            'results.frames.csv',
            '`GMXMMPBSA_radii.json`',
            'diagnostic zip bundle',
        )
        for text in required:
            with self.subTest(text=text):
                self.assertIn(text, guide)

    def test_changelog_links_to_170_migration_guide(self):
        changelog = (ROOT / 'docs/changelog.md').read_text()
        self.assertIn('compatibility.md#migrating-from-165-to-170', changelog)


class FeatureOverviewDocumentationTest(unittest.TestCase):
    def test_overviews_list_all_console_applications(self):
        for filename in ('docs/summary.md', 'docs/getting-started.md'):
            document = (ROOT / filename).read_text()
            with self.subTest(filename=filename):
                self.assertIn('gmx_MMPBSA', document)
                self.assertIn('amber_MMPBSA', document)
                self.assertIn('gmx_MMPBSA_ana', document)
                self.assertIn('gmx_MMPBSA_test', document)
                self.assertIn('four applications', document)

    def test_overviews_describe_composite_mutation_semantics(self):
        overview = (ROOT / 'docs/summary.md').read_text()
        alanine = (ROOT / 'docs/examples/Alanine_scanning/README.md').read_text()
        for document in (overview, alanine):
            normalized = ' '.join(document.split())
            with self.subTest(document=document[:30]):
                self.assertIn('one residue or several residues together in one composite mutant', normalized)
                self.assertIn('does not independently scan each selected residue', normalized)
                self.assertNotIn('Only one residue can be mutated per calculation', normalized)

    def test_changelog_matches_native_amber_st_and_mt_scope(self):
        changelog = (ROOT / 'docs/changelog.md').read_text()
        self.assertIn('supported ST and MT paths', changelog)
        self.assertIn('IE/C2 with MT is experimental', changelog)
        self.assertNotIn('ligand MT unsupported', changelog)


class CommandReferenceDocumentationTest(unittest.TestCase):
    def test_gmx_command_reference_matches_current_help_contract(self):
        document = ' '.join((ROOT / 'docs/gmx_MMPBSA_command-line.md').read_text().split())
        for text in (
                'gbnsr6',
                'Allowed formats: *.tpr (recommended), *.pdb',
                'GMXMMPBSA/data/xvv_files/tip3p.xvv',
                '--no-error-bundle',
                '2026.0'):
            with self.subTest(text=text):
                self.assertIn(text, document)
        self.assertNotIn('*.gro (default: None)', document)

    def test_amber_command_reference_covers_shared_and_native_options(self):
        document = (ROOT / 'docs/amber_MMPBSA.md').read_text()
        for text in (
                'amber_MMPBSA -h',
                '--create_input',
                '-nogui',
                '--no-error-bundle',
                '--rewrite-output',
                '-cp',
                '-cm',
                '-ct',
                '-rp/-rm/-rt',
                '-lp/-lm/-lt'):
            with self.subTest(text=text):
                self.assertIn(text, document)


class InputReferenceDocumentationTest(unittest.TestCase):
    def test_documented_defaults_and_syntax_match_parser(self):
        document = (ROOT / 'docs/input_file.md').read_text()
        for text in (
                '`treeCoulomb` (Default = 0)',
                '`idecomp` (Default = 2)',
                '`dec_verbose` (Default = 1)',
                '`nmstartframe`[^2] (Default = 1)',
                '`solvcut` (Default = -1)',
                '`ndiis_attempts` (Default = None)',
                'The accepted levels are 0, 1, and 2.',
                '`#` starts an inline comment inside a'):
            with self.subTest(text=text):
                self.assertIn(text, document)

    def test_noasympcorr_documentation_matches_parser_semantics(self):
        document = (ROOT / 'docs/input_file.md').read_text()
        variable = input_file.namelists['rism'].variables['noasympcorr']

        self.assertEqual(variable.value, 1)
        self.assertIn('Disable asymptotic corr.', variable.description)
        self.assertIn('Disable long-range asymptotic corrections', document)
        self.assertIn('Long-range asymptotics are still used to', document)
        self.assertNotIn('1: Use the long-range corrections', document)


class RismAndStatisticsDocumentationTest(unittest.TestCase):
    def test_rism_docs_match_sander_execution_and_current_artifacts(self):
        input_doc = (ROOT / 'docs/input_file.md').read_text()
        output_doc = (ROOT / 'docs/output.md').read_text()

        self.assertIn('AmberTools `sander` backend', input_doc)
        self.assertIn('historical pre-v1.5.2 backend reference', input_doc)
        self.assertIn('`_GMXMMPBSA_rism.mdin`', output_doc)
        self.assertIn('`_GMXMMPBSA_complex_rism.mdout.#`', output_doc)
        self.assertIn('frame distribution by', output_doc)
        self.assertNotIn('`_GMXMMPBSA_complex_rism.out.#`', output_doc)

    def test_block_statistics_docs_state_algorithm_and_limitations(self):
        document = (ROOT / 'docs/output.md').read_text()
        for text in (
                'candidate block sizes are',
                'at least eight',
                'trailing frames that do not fill a complete block are excluded',
                '`ddof=1`',
                'legacy population frame SEM',
                'does not establish convergence',
                'force-field, solvent-model, or other model error',
                'Frame `ddof=0` vs block `ddof=1`',
                'Average** column is unaffected'):
            with self.subTest(text=text):
                self.assertIn(text, document)


class ApiAndNotebookDocumentationTest(unittest.TestCase):
    def test_api_docs_and_example_include_both_block_summary_rows(self):
        api_doc = (ROOT / 'docs/api.md').read_text()
        example_doc = (ROOT / 'examples/API/README.md').read_text()
        script = (ROOT / 'examples/API/extract_api_data.py').read_text()
        for document in (api_doc, example_doc, script):
            with self.subTest(document=document[:30]):
                self.assertIn('Block SD', document)
                self.assertIn('Block SEM', document)

    def test_local_notebook_is_safe_by_default_and_filters_numeric_frames(self):
        notebook = (ROOT / 'notebooks/gmx_MMPBSA_Local.ipynb').read_text()
        self.assertIn('RUN_LOCAL_CALC = False', notebook)
        self.assertIn('physical temporary working directory', notebook)
        self.assertIn('numeric_frame_mask', notebook)
        self.assertIn('Energy (kcal/mol)', notebook)
        self.assertIn('UNCERTAINTY = \\"Block SEM\\"', notebook)
        self.assertNotIn('RUN_LOCAL_CALC = True', notebook)


class TesterAndScientificDocumentationTest(unittest.TestCase):
    def test_tester_docs_isolate_local_examples_and_match_current_help(self):
        document = (ROOT / 'docs/examples/gmx_MMPBSA_test.md').read_text()
        for text in (
                'cp -a ./examples/.',
                '--examples-dir "$TMP_EXAMPLES"',
                'physical temporary tree',
                '`-r/--reuse`',
                '8    x | 10  Metalloprotein-ligand',
                '15   . | 10  Interaction Entropy approximation',
                '17   x | 10  Entropy calculation using Normal Mode approximation'):
            with self.subTest(text=text):
                self.assertIn(text, document)
        self.assertNotIn('--examples-dir ./examples -t 2', document)

    def test_active_multicomponent_example_is_visible(self):
        examples_index = (ROOT / 'examples/README.md').read_text()
        nav = (ROOT / 'mkdocs.yml').read_text()
        self.assertIn('[Multicomponent system (Comp_receptor)](Comp_receptor/README.md)', examples_index)
        self.assertIn('- Multicomponent system: examples/Comp_receptor/README.md', nav)

    def test_pb_rewrite_and_entropy_docs_preserve_scientific_limits(self):
        qa = (ROOT / 'docs/Q&A/calculations.md').read_text()
        introduction = (ROOT / 'docs/introduction.md').read_text()
        output = (ROOT / 'docs/output.md').read_text()

        for text in (
                'legacy post-processing workaround',
                'does not rerun PB',
                'recalculated `inp=1`',
                'do not edit the original result in place'):
            with self.subTest(text=text):
                self.assertIn(text, qa)
        for text in (
                'mass-weighted Hessian',
                'coordinate covariance matrix',
                'not universally superior to NMODE',
                'not generally sufficient by itself for relative affinity claims'):
            with self.subTest(text=text):
                self.assertIn(text, introduction)
        for text in ('`1-4 VDW`', '`1-4 EEL`', '`UB`', '`IMP`', '`CMAP`', '`ESCF`',
                     'not automatically an interaction energy'):
            with self.subTest(text=text):
                self.assertIn(text, output)


class ReleasePublicationDocumentationTest(unittest.TestCase):
    def test_output_catalog_distinguishes_legacy_and_portable_artifacts(self):
        output = (ROOT / 'docs/output.md').read_text()
        input_doc = (ROOT / 'docs/input_file.md').read_text()
        for text in (
                'Legacy output illustration',
                'historical v1.4.3-format illustration',
                '`COMPACT_MMXSA_RESULTS.mmxsa`',
                '`GMXMMPBSA_radii.json`',
                '`gmx_MMPBSA_error_bundle_*.zip`',
                'Legacy QH artifacts'):
            with self.subTest(text=text):
                self.assertIn(text, output)
        self.assertIn('_COMPACT_MMXSA_RESULTS.mmxsa_', input_doc)
        self.assertNotIn('_COMPACT_gmx_MMPBSA_RESULTS.mmxsa_', input_doc)

    def test_proposed_support_claims_are_bounded(self):
        readme = (ROOT / 'README.md').read_text()
        getting_started = (ROOT / 'docs/getting-started.md').read_text()
        versus = (ROOT / 'docs/versus.md').read_text()
        howworks = (ROOT / 'docs/howworks.md').read_text()
        self.assertIn('GROMACS `>=2022,<2027`', readme)
        self.assertIn('does not mean that every system', getting_started)
        self.assertIn('legacy QH reader', versus)
        self.assertIn('many supported [MMPBSA.py]', howworks)
        self.assertNotIn('works with all GROMACS versions', readme.lower())

    def test_canonical_example_links_target_repository_docs(self):
        root_readme = (ROOT / 'README.md').read_text()
        protein = (ROOT / 'examples/Protein_protein/README.md').read_text()
        nested = (ROOT / 'examples/Protein_ligand/ST/README.md').read_text()
        self.assertIn('(docs/cite_us.md)', root_readme)
        self.assertIn('../../docs/input_file.md', protein)
        self.assertIn('../../../docs/input_file.md', nested)
        self.assertIn('../../docs/examples/gmx_MMPBSA_test.md', protein)
        self.assertIn('../../../docs/examples/gmx_MMPBSA_test.md', nested)

    def test_publication_inventory_has_explicit_decisions(self):
        config = (ROOT / 'mkdocs.yml').read_text()
        explicit_waters = (ROOT / 'examples/Explicit_receptor_waters/README.md').read_text()
        changelog = (ROOT / 'docs/changelog.md').read_text()
        getting_started = (ROOT / 'docs/getting-started.md').read_text()
        for text in (
                'application_note_advances_since_1_4_3.md',
                'spikes/**',
                'examples/COVID-19_related_proteins/**',
                'examples/Protein_DNA_RNA_Ion_ligand/**',
                'Support: support.md'):
            with self.subTest(text=text):
                self.assertIn(text, config)
        app_note = (ROOT / 'docs/application_note_advances_since_1_4_3.md').read_text()
        self.assertIn('frozen through v1.6.5', app_note)
        self.assertIn('1.6.5→1.7.0 migration guide', app_note)
        self.assertNotIn('summarizes advances through v1.6.5', app_note)
        self.assertIn('versioned example archive', explicit_waters)
        self.assertNotIn('downgit.github.io/#/home?url=https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tree/master/examples/Explicit_receptor_waters', explicit_waters)
        self.assertIn('examples/QM_MMGBSA/README.md', changelog)
        self.assertNotIn('sourcery.ai/pro', (ROOT / 'README.md').read_text())
        self.assertNotIn('sourcery.ai/pro', getting_started)


class DocumentationInputExamplesTest(unittest.TestCase):
    def _fenced_block(self, document, predicate):
        blocks = re.findall(r'```[^\n]*\n(.*?)```', document, flags=re.DOTALL)
        block = next((candidate for candidate in blocks if predicate(candidate)), None)
        self.assertIsNotNone(block)
        return block

    def _parse_fenced_block(self, document, predicate):
        block = self._fenced_block(document, predicate)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'mmpbsa.in'
            path.write_text(block)
            parser = copy.deepcopy(input_file)
            for namelist in parser.namelists.values():
                namelist.open = False
            return parser.Parse(path)

    def test_rism_documentation_example_uses_current_correction_controls(self):
        document = (ROOT / 'docs/input_file.md').read_text()
        block = self._fenced_block(document, lambda candidate: '&rism' in candidate)
        parsed = self._parse_fenced_block(document, lambda candidate: '&rism' in candidate)

        self.assertEqual(parsed['rism']['polardecomp'], 1)
        self.assertEqual(parsed['rism']['gfcorrection'], 1)
        self.assertNotIn('thermo', block)
        self.assertNotIn('rgbmax', document)

    def test_citation_documentation_example_is_a_complete_namelist(self):
        document = (ROOT / 'docs/cite_us.md').read_text()
        parsed = self._parse_fenced_block(document, lambda block: '&decomp' in block)

        self.assertEqual(parsed['decomp']['idecomp'], 2)
        self.assertEqual(parsed['decomp']['dec_verbose'], 3)
        self.assertEqual(parsed['decomp']['print_res'], 'within 5')


if __name__ == '__main__':
    unittest.main()
