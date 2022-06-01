from typing import Union
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import datad
from pathlib import Path
import json
from GMXMMPBSA.analyzer.style import *


def flatten_dict(d, parent_key=None):
    if parent_key is None:
        parent_key = []
    items = []
    for k, v in d.items():
        new_key = tuple(list(parent_key) + [k if parent_key else k])
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key).items())
        elif k not in ['value', 'action_type', 'default'] or v is None:
            continue
        else:
            items.append((new_key, v))
    return dict(items)


def value2deafult(d):
    """
    Iterate recursively the settings dict and make default = value
    """
    if isinstance(d, dict):
        if 'default' in d and 'value' in d:
            d['default'] = d['value']
        if 'children' in d:
            value2deafult(d['children'])
        else:
            for k, v in d.items():
                value2deafult(v)


tooltip1 = '''
            <html lang="en">

<style>
    table, th {
        border-bottom: 2px solid black;
        /*border-top: 2px solid black;*/
    }
    img {
        width: 300px;
        height: 25px;
    }
</style>
<body>
<h3>Palettes</h3>
<p>Tooltip</p>''' + f'''
<table>
     <tr>
        <th colspan="5">Cartoon Colors</th>
    </tr>
    <tr>
        <td>ArmyRose_7</td>
        <td><img src='{ArmyRose_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>ArmyRose_5</td>
        <td><img src='{ArmyRose_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>Geyser_7</td>
        <td><img src='{Geyser_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>Geyser_5</td>
        <td><img src='{Geyser_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>TealRose_7</td>
        <td><img src='{TealRose_7}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>TealRose_5</td>
        <td><img src='{TealRose_5}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>Tropic_7</td>
        <td><img src='{Tropic_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>Tropic_5</td>
        <td><img src='{Tropic_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <th colspan="5">CMOcean Colors</th>
    </tr>
    <tr>
        <td>Balance_7</td>
        <td><img src='{Balance_7}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>Balance_5</td>
        <td><img src='{Balance_5}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>Curl_7</td>
        <td><img src='{Curl_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>Curl_5</td>
        <td><img src='{Curl_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>Delta_7</td>
        <td><img src='{Delta_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>Delta_5</td>
        <td><img src='{Delta_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <th colspan="5">Colorbrew Colors</th>
    </tr>
    <tr>
        <td>BrBG_7</td>
        <td><img src='{BrBG_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>BrBG_5</td>
        <td><img src='{BrBG_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>PiYG_7</td>
        <td><img src='{PiYG_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>PiYG_5</td>
        <td><img src='{PiYG_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>PRGn_7</td>
        <td><img src='{PRGn_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>PRGn_5</td>
        <td><img src='{PRGn_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>PuOr_7</td>
        <td><img src='{PuOr_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>PuOr_5</td>
        <td><img src='{PuOr_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>RdBu_7</td>
        <td><img src='{RdBu_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>RdBu_5</td>
        <td><img src='{RdBu_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>RdGy_7</td>
        <td><img src='{RdGy_7}' width="200" height="16" alt=""></td>
        <td></td>
        <td>RdGy_5</td>
        <td><img src='{RdGy_5}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>RdYlBu *</td>
        <td><img src='{RdYlBu}' width="200" height="16" alt=""></td>
        <td></td>
        <td>RdYlGn *</td>
        <td><img src='{RdYlGn}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>Spectral *</td>
        <td><img src='{Spectral}' width="200" height="16" alt=""></td>
        <td></td>
        <td>coolwarm *</td>
        <td><img src='{coolwarm}' width="200" height="16" alt=""></td>
    </tr>
</table>

<p>* No compatible with PyMOL visualization</p>
</body>
</html>'''
tooltip2 = '''
            <html lang="en">

<style>
    table, th {
        border-bottom: 2px solid black;
        /*border-top: 2px solid black;*/
    }
    img {
        width: 300px;
        height: 25px;
    }
</style>
<body>
<h3>Palettes</h3>
<p>Tooltip</p>''' + f'''
<table>
     <tr>
        <th colspan="5">Color cycle - 20 -</th>
    </tr>
    <tr>
        <td>husl</td>
        <td><img src='{husl}' width="200" 
        height="16" alt=""></td>
        <td></td>
        <td>hls</td>
        <td><img src='{hls}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>tab20</td>
        <td><img src='{tab20}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>tab20b</td>
        <td><img src='{tab20b}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>tab20c</td>
        <td><img src='{tab20c}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <th colspan="5">Color cycle - 12 -</th>
    </tr>
    <tr>
        <td>Paired</td>
        <td><img src='{Paired}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>Set3</td>
        <td><img src='{Set3}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <th colspan="5">Color cycle - 10 -</th>
    </tr>
    <tr>
        <td>deep</td>
        <td><img src='{deep}' width="200" height="16" alt=""></td>
        <td></td>
        <td>muted</td>
        <td><img src='{muted}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <td>bright</td>
        <td><img src='{bright}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>pastel</td>
        <td><img src='{pastel}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>dark</td>
        <td><img src='{dark}' width="200" height="16" alt=""></td>
        <td></td>
        <td>colorblind</td>
        <td><img src='{colorblind}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>tab10</td>
        <td><img src='{tab10}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <th colspan="5">Color cycle - 9 -</th>
    </tr>
    <tr>
        <td>Pastel1</td>
        <td><img src='{Pastel1}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>Set1</td>
        <td><img src='{Set1}' width="200" height="16" alt=""></td>
    </tr>
    <tr>
        <th colspan="5">Color cycle - 8 -</th>
    </tr>
    <tr>
        <td>Pastel2</td>
        <td><img src='{Pastel2}' width="200" height="16" 
        alt=""></td>
        <td></td>
        <td>Accent</td>
        <td><img src='{Accent}' width="200" height="16" 
        alt=""></td>
    </tr>
    <tr>
        <td>Dark2</td>
        <td><img src='{Dark2}' width="200" height="16" alt=""></td>
        <td></td>
        <td>Set2</td>
        <td><img src='{Set2}' width="200" height="16" alt=""></td>
    </tr>
</table>
<p>Please</p>
</body>
</html>'''


class ChartSettings(dict):
    D = 1  # Drawable
    U = 2  # Updatable
    R = 3  # re-plot
    default = {
        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'parameters', 'value': None, 'default': None,
        'children': {
            'General': {
                'type': 'group', 'enabled': True, 'expanded': True, 'name': 'General', 'value': None, 'default':
                    None,
                'children': {
                    'theme': {
                        'type': 'list', 'enabled': True, 'expanded': True, 'name': 'theme', 'value': 'darkgrid',
                        'values': ['darkgrid', '---0', 'whitegrid', 'dark', 'white', 'ticks'], 'default': 'darkgrid',
                        'action_type': R, 'tip': 'changes visual theme for all plots'},
                    'toolbar': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'toolbar', 'value': False,
                                'default': False, 'action_type': U, 'tip': 'shows toolbar in all graphs'},
                    'figure-format': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'figure-format', 'value': None,
                        'default': None,
                        'children': {
                            'dpi-plot': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'dpi-plot',
                                         'value': 100, 'step': 10, 'accelerated': True, 'limits': (50, 300),
                                         'default': 100, 'action_type': R, 'tip': 'set dpi for all plots'},
                            'dpi-save': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'dpi-save',
                                         'value': 300, 'step': 10, 'accelerated': True, 'limits': (50, 1200),
                                         'default': 300, 'action_type': U, 'tip': 'dpi used for saved figures'},
                            'save-format': {
                                'type': 'list', 'enabled': True, 'expanded': True, 'name': 'save-format',
                                'value': 'svg', 'values': ['svg', 'png', 'tiff', 'pdf', 'eps', 'jpg', 'jpeg', 'pgf',
                                                           'ps', 'raw', 'rgba', 'svgz', 'tif'],
                                'default': 'svg', 'action_type': U},
                        }}}},
            'Line Plot': {
                'type': 'group', 'enabled': True, 'expanded': False, 'name': 'Line Plot', 'value': None, 'default':
                    None,
                'children': {
                    'line-width': {'type': 'float', 'enabled': True, 'expanded': True, 'name': 'line-width',
                                   'value': 0.7, 'step': 0.1, 'limits': (0.1, 1.5), 'accelerated': True,
                                   'default': 0.7, 'action_type': R},
                    'line-color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'line-color',
                                   'value': [0, 0, 0, 255], 'default': [0, 0, 0, 255], 'action_type': R},

                    'Rolling average': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'Rolling average', 'value': None,
                        'default': None,
                        'children': {
                            'show': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'show', 'value': True,
                                     'default': True, 'action_type': R},
                            'color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'color',
                                      'value': [255, 0, 0, 255], 'default': [255, 0, 0, 255], 'action_type': R},
                            'style': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'style',
                                      'value': 'dashed', 'values': ['dashed', 'solid', 'dashdot', 'dotted'],
                                      'default': 'dashed', 'action_type': R},
                            'width': {'type': 'float', 'enabled': True, 'expanded': True, 'name': 'width',
                                      'value': 0.8, 'step': 0.1, 'default': 0.8, 'action_type': R},
                            'window': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'window',
                                       'value': 50, 'step': 1, 'limits': (2, 100000), 'accelerated': True,
                                       'default': 50, 'action_type': R},
                            'first_obs': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'first_obs',
                                          'value': 'window', 'values': ['window', 'start'], 'default': 'window',
                                          'action_type': R, 'tip': 'Defines the start value of the moving average. '
                                                                   'Window: set start according to'
                                                                   'the defined window size. Start: '
                                                                   'set start at the first value.'}}},
                    'fontsize': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'fontsize', 'value': None,
                        'default': None,
                        'children': {
                            'x-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'y-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'x-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'y-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'title': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'title',
                                      'value': 14, 'default': 14, 'action_type': D},
                            'suptitle': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'suptitle',
                                         'value': 12, 'default': 12, 'action_type': D},
                            'legend': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'legend',
                                       'value': 9, 'default': 9, 'action_type': D}}},
                    'axes': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'axes', 'value': None,
                        'default': None,
                        'children': {
                            'num-xticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'num-xticks',
                                           'value': 10, 'default': 10, 'action_type': R},
                            'num-yticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'num-yticks',
                                           'value': 10, 'default': 10, 'action_type': R},
                            'x-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-rotation',
                                           'value': 0, 'step': 1, 'limits': (-90, 90), 'accelerated': True,
                                           'default': 0, 'action_type': D},
                            'y-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-rotation',
                                           'value': 0, 'step': 1, 'limits': (-90, 90), 'accelerated': True,
                                           'default': 0, 'action_type': D}}},
                    'figure': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'figure', 'value': None,
                        'default': None,
                        'children': {
                            'width': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'width', 'value': 8,
                                      'default': 8, 'action_type': R},
                            'height': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'height',
                                       'value': 4, 'default': 4, 'action_type': R}}}}},
            'Bar Plot': {
                'type': 'group', 'enabled': True, 'expanded': False, 'name': 'Bar Plot', 'value': None, 'default': None,
                'children': {
                    'use-palette': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'use-palette',
                                    'value': True, 'default': True, 'action_type': R},
                    'palette': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'palette',
                                'value': 'husl', 'values': ['husl', '---0', 'hls', 'deep', 'muted', 'bright', 'pastel',
                                                            'dark', 'colorblind', '---1',
                                                            'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2', 'Set1',
                                                            'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c'],
                                'default': 'husl', 'action_type': R, 'tip': tooltip2},
                    'color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'color',
                              'value': [44, 105, 176, 255], 'default': [44, 105, 176, 255], 'action_type': R},
                    'subplot-components': {'type': 'bool', 'enabled': True, 'expanded': True,
                                           'name': 'subplot-components', 'value': True, 'default': True,
                                           'action_type': R, 'tip': 'Groups energy components in bar plot'},
                    'scale-yaxis': {'type': 'bool', 'enabled': True, 'expanded': True,
                                    'name': 'scale-yaxis', 'value': False, 'default': False, 'action_type': R,
                                    'tip': 'Uses exponential values in y-axis'},
                    'remove-molid': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'remove-molid',
                                     'value': True, 'default': True, 'action_type': R},
                    'IE/C2 Entropy': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'IE/C2 Entropy', 'value': None,
                        'default': None,
                        'children': {
                            'ie-color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'ie-color',
                                         'value': [0, 0, 255, 255], 'default': [0, 0, 255, 255], 'action_type': R},
                            'sigma-color': {
                                'type': 'group', 'enabled': True, 'expanded': True, 'name': 'sigma-color',
                                'value': None, 'default': None,
                                'children': {
                                    'reliable': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'reliable',
                                                 'value': [0, 255, 0, 255], 'default': [0, 255, 0, 255],
                                                 'action_type': R, 'tip': 'Color if sigma < 3.6kcal/mol'},
                                    'non-reliable': {'type': 'color', 'enabled': True, 'expanded': True,
                                                     'name': 'non-reliable', 'value': [255, 0, 0, 255],
                                                     'default': [255, 0, 0, 255], 'action_type': R,
                                                     'tip': 'Color if sigma > 3.6kcal/mol'}}}}},
                    'error-line': {'type': 'group', 'enabled': True, 'expanded': True, 'name': 'error-line',
                                   'value': None, 'default': None,
                                   'children': {
                                       'width': {'type': 'float', 'enabled': True, 'expanded': True, 'name': 'width',
                                                 'value': 0.7, 'step': 0.1, 'limits': (0.1, 1.5),
                                                 'accelerated': True, 'default': 0.7, 'action_type': R},
                                       'color': {'type': 'color', 'enabled': True, 'expanded': True,
                                                 'name': 'color', 'value': [0, 0, 0, 255],
                                                 'default': [0, 0, 0, 255], 'action_type': R},
                                       'cap-size': {'type': 'int', 'enabled': True, 'expanded': True,
                                                    'name': 'cap-size', 'value': 0, 'step': 1, 'limits': (0, 50),
                                                    'accelerated': True, 'default': 0, 'action_type': R},
                                       'representation': {'type': 'list', 'enabled': True, 'expanded': True,
                                                          'name': 'representation', 'value': 'SD',
                                                          'values': ['SD', 'SEM'], 'default': 'SD', 'action_type': R,
                                                          'tip': "The metric representing the error line. SD, the "
                                                                 "Standard Deviation, describes the data distribution "
                                                                 "around the mean, while SEM, the Standard Error of "
                                                                 "the Mean, describes how representative the mean is "
                                                                 "of the population."},


                                   }},
                    'bar-label': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'bar-label', 'value': None,
                        'default': None,
                        'children': {
                            'show': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'show', 'value': False,
                                     'default': False, 'action_type': R},
                            'fontsize': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'fontsize',
                                         'value': 8, 'limits': (2, 20), 'default': 8, 'action_type': D},
                            'label_type': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'label_type',
                                           'value': 'edge', 'values': ['edge', 'center'], 'default': 'edge',
                                           'action_type': R},
                            'padding': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'padding',
                                        'value': 5, 'limits': (0, 50), 'default': 5, 'action_type': R}
                        }},
                    'fontsize': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'fontsize', 'value': None,
                        'default': None,
                        'children': {
                            'x-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'y-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'x-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'y-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'title': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'title',
                                      'value': 14, 'default': 14, 'action_type': D},
                            'suptitle': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'suptitle',
                                         'value': 12, 'default': 12, 'action_type': D}}},
                    'axes': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'axes', 'value': None,
                        'default': None,
                        'children': {
                            'num-yticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'num-yticks',
                                           'value': 10, 'default': 10, 'action_type': D},
                            'x-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-rotation',
                                           'value': 45, 'step': 1, 'limits': (-90, 90),
                                           'accelerated': True, 'default': 45, 'action_type': D},
                            'y-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-rotation',
                                           'value': 0, 'step': 1, 'limits': (-90, 90),
                                           'accelerated': True, 'default': 0, 'action_type': D},
                            'y-inverted': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'y-inverted',
                                           'value': True, 'default': True, 'action_type': R}
                        }},
                    'figure': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'figure', 'value': None,
                        'default': None,
                        'children': {
                            'width': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'width',
                                      'value': 5, 'default': 5, 'action_type': R},
                            'height': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'height',
                                       'value': 4, 'default': 4, 'action_type': R}}}}},
            'Heatmap Plot': {
                'type': 'group', 'enabled': True, 'expanded': False, 'name': 'Heatmap Plot', 'value': None,
                'default': None,
                'children': {
                    'highlight-components': {'type': 'bool', 'enabled': True, 'expanded': True,
                                             'name': 'highlight-components', 'value': True, 'default': True,
                                             'action_type': R, 'tip': 'Highlights receptor/ligand residues'},
                    'legend': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'legend', 'value': True,
                               'default': True, 'action_type': R},
                    'remove-molid': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'remove-molid',
                                     'value': True, 'default': True, 'action_type': R,
                                     'tip': 'Removes R/L labels from residue names'},
                    'receptor-color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'receptor-color',
                                       'value': [255, 85, 127, 255], 'default': [255, 85, 127, 255], 'action_type': R},
                    'ligand-color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'ligand-color',
                                     'value': [0, 255, 0, 255], 'default': [0, 255, 0, 255], 'action_type': R},

                    'y-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-rotation',
                                   'value': 0, 'step': 1, 'limits': (-90, 90), 'accelerated': True, 'default': 0,
                                   'action_type': D},
                    'Per-wise': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'Per-wise', 'value': None,
                        'default': None,
                        'children': {
                            'palette': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'palette',
                                        'value': 'seismic', 'values': ['seismic', '---1',
                                                                       'Balance_7', 'Balance_5', 'Curl_7', 'Curl_5',
                                                                       'Delta_7', 'Delta_5', '---2',
                                                                       'ArmyRose_7', 'ArmyRose_5', 'Geyser_7',
                                                                       'Geyser_5', 'Temps_7', 'Temps_5', 'TealRose_7',
                                                                       'TealRose_5', 'Tropic_7', 'Tropic_5', '---3',
                                                                       'BrBG_7', 'BrBG_5', 'PRGn_7', 'PRGn_5', 'PiYG_7',
                                                                       'PiYG_5', 'PuOr_7', 'PuOr_5', 'RdBu_7', 'RdBu_5',
                                                                       'RdGy_7', 'RdGy_5', '---4', '---5',
                                                                       'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm'
                                                                       ], 'default': 'seismic', 'action_type': R,
                                        'tip': tooltip1},
                            'annotation': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'annotation',
                                           'value': False, 'default': False, 'action_type': R},
                            'x-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-rotation',
                                           'value': 90, 'step': 1, 'limits': (-90, 90), 'accelerated': True,
                                           'default': 90, 'action_type': D}}},
                    'Per-residue': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'Per-residue', 'value': None,
                        'default': None,
                        'children': {
                            'palette': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'palette',
                                        'value': 'seismic', 'values': ['seismic', '---1',
                                                                       'Balance_7', 'Balance_5', 'Curl_7', 'Curl_5',
                                                                       'Delta_7', 'Delta_5', '---2',
                                                                       'ArmyRose_7', 'ArmyRose_5', 'Geyser_7',
                                                                       'Geyser_5', 'Temps_7', 'Temps_5', 'TealRose_7',
                                                                       'TealRose_5', 'Tropic_7', 'Tropic_5', '---3',
                                                                       'BrBG_7', 'BrBG_5', 'PRGn_7', 'PRGn_5', 'PiYG_7',
                                                                       'PiYG_5', 'PuOr_7', 'PuOr_5', 'RdBu_7', 'RdBu_5',
                                                                       'RdGy_7', 'RdGy_5', '---4', '---5',
                                                                       'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm'
                                                                       ], 'default': 'seismic', 'action_type': R,
                                        'tip': tooltip1},
                            'num-xticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'num-xticks',
                                           'value': 10, 'default': 10, 'action_type': R},
                            'x-rotation': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-rotation',
                                           'value': 0, 'step': 1, 'limits': (-90, 90), 'accelerated': True,
                                           'default': 0, 'action_type': D}}},
                    'fontsize': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'fontsize', 'value': None,
                        'default': None,
                        'children': {
                            'x-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'y-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'x-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'y-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'title': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'title',
                                      'value': 14, 'default': 14, 'action_type': D},
                            'suptitle': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'suptitle',
                                         'value': 12, 'default': 12, 'action_type': D},
                            'legend': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'legend',
                                       'value': 9, 'default': 9, 'action_type': R},
                            'colorbar-ticks': {'type': 'int', 'enabled': True, 'expanded': True,
                                               'name': 'colorbar-ticks', 'value': 9, 'default': 9, 'action_type': R},
                            'colorbar-label': {'type': 'int', 'enabled': True, 'expanded': True,
                                               'name': 'colorbar-label', 'value': 11, 'default': 11, 'action_type': R},
                            'annotation': {'type': 'int', 'enabled': True, 'expanded': True,
                                           'name': 'annotation', 'value': 8, 'default': 8, 'action_type': R}
                        }},
                    'figure': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'figure', 'value': None,
                        'default': None,
                        'children': {
                            'width-per-residue': {'type': 'int', 'enabled': True, 'expanded': True,
                                                  'name': 'width-per-residue', 'value': 10, 'default': 10,
                                                  'action_type': R},
                            'width-per-wise': {'type': 'int', 'enabled': True, 'expanded': True,
                                               'name': 'width-per-wise', 'value': 8, 'default': 8, 'action_type': R},
                            'height': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'height',
                                       'value': 7, 'default': 7, 'action_type': R}}}}},
            'Visualization': {
                'type': 'group', 'enabled': True, 'expanded': False, 'name': 'Visualization', 'value': None,
                'default': None,
                'children': {
                    'palette': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'palette',
                                'value': 'auto', 'values': ['auto', '---0', 'seismic', '---1',
                                                            'Balance_7', 'Balance_5', 'Curl_7', 'Curl_5',
                                                            'Delta_7', 'Delta_5', '---2',
                                                            'ArmyRose_7', 'ArmyRose_5', 'Geyser_7',
                                                            'Geyser_5', 'Temps_7', 'Temps_5', 'TealRose_7',
                                                            'TealRose_5', 'Tropic_7', 'Tropic_5', '---3',
                                                            'BrBG_7', 'BrBG_5', 'PRGn_7', 'PRGn_5', 'PiYG_7',
                                                            'PiYG_5', 'PuOr_7', 'PuOr_5', 'RdBu_7', 'RdBu_5',
                                                            'RdGy_7', 'RdGy_5'
                                                            ], 'default': 'auto', 'action_type': R, 'tip': tooltip1},
                    'cartoon_oval_length': {'type': 'float', 'enabled': True, 'expanded': True,

                                            'name': 'cartoon_oval_length', 'value': 1.0, 'step': 0.1, 'default': 1.0,
                                            'action_type': R},
                    'cartoon_rect_length': {'type': 'float', 'enabled': True, 'expanded': True,
                                            'name': 'cartoon_rect_length', 'value': 1.2, 'step': 0.1, 'default': 1.2,
                                            'action_type': R},
                    'cartoon_rect_width': {'type': 'float', 'enabled': True, 'expanded': True,
                                           'name': 'cartoon_rect_width', 'value': 0.3, 'step': 0.1, 'default': 0.3,
                                           'action_type': R},
                    'cartoon_side_chain_helper': {'type': 'bool', 'enabled': True, 'expanded': True,
                                                  'name': 'cartoon_side_chain_helper', 'value': True, 'default': True,
                                                  'action_type': R},
                    'light_count': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'light_count',
                                    'value': 1, 'step': 1, 'accelerated': False, 'limits': (1, 10),
                                    'default': 1, 'action_type': R, 'tip': 'Defines the number of light sources'},
                    'background': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'background',
                                   'value': 'gray50', 'values': ['gray50', '---0', 'black',
                                                                 'gray10', 'gray20', 'gray30', 'gray40',
                                                                 'gray60', 'gray70', 'gray80', 'gray90', 'white'],
                                   'default': 'gray50', 'action_type': R},
                    'representation': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'representation',
                                       'value': 'sticks', 'values': ['sticks', '---0', 'lines', '---1', 'spheres',
                                                                     'surface', 'mesh', 'dots', 'lines+dots',
                                                                     'sticks+dots'],
                                       'default': 'sticks', 'action_type': R},
                }}
        }}

    def __init__(self, custom: Union[Path, str] = None):
        super(ChartSettings, self).__init__()

        self.config_folder = Path('~').expanduser().absolute().joinpath('.config', 'gmx_MMPBSA')
        self.filename = self.config_folder.joinpath('settings.json')

        if isinstance(custom, Path):
            with open(custom.joinpath('settings.json')) as read_file:
                config = json.load(read_file)
                self.update(config)
        elif isinstance(custom, str):
            with open(self.filename) as read_file:
                config = json.load(read_file)
                self.update(config)
        else:
            self.update(self.default)

        self.changes = dict(
            line_action=0,
            line_ie_action=0,
            bar_action=0,
            heatmap_action=0,
            visualization_action=0,
            figure=0
        )

    def config_created(self):
        return self.filename.exists()

    def return_default(self):
        self.update(self.default)

    def write_system_config(self, syspath: Path = None):
        """

        @param syspath: Current system path
        @return:
        """
        if syspath is None:
            self.config_folder.mkdir(exist_ok=True)
            with open(self.filename, "w") as write_file:
                json.dump(self, write_file, indent=4)
        else:
            self.set_as_default()
            filename = syspath.joinpath('settings.json')
            with open(filename, "w") as write_file:
                json.dump(self, write_file, indent=4)

    def _get_change_type(self, rdict_value, var):
        if rdict_value == self.D:
            self.changes[f'{var}_redraw'] = True
        elif rdict_value == self.R:
            self.changes[f'{var}_replot'] = True
        else:
            self.changes[f'{var}_update'] = True

    def is_changed(self, osett) -> bool:
        flatten = flatten_dict(self)
        f_osett = flatten_dict(osett)
        return flatten != f_osett

    def is_default_changed(self):
        flatten = flatten_dict(self)
        for x, v in flatten.items():
            if x[-1] != 'value':
                continue
            ks = list(x[:-1])
            if v != flatten[tuple(ks + ['default'])]:
                return True

    def get_changes(self, osett) -> None:

        flatten = flatten_dict(self)
        f_osett = flatten_dict(osett)
        if flatten == f_osett:
            return
        for k, v in flatten.items():
            if k[-1] == 'value' and v != f_osett[k]:
                at = tuple(list(k[:-1]) + ['action_type'])
                if k[1] == 'Bar Plot':
                    self.changes['bar_action'] = flatten[at]
                if k[3] == 'IE/C2 Entropy':
                    self.changes['line_ie_action'] = flatten[at]
                elif k[1] == 'Line Plot':
                    self.changes['line_action'] = flatten[at]
                elif k[1] == 'Heatmap Plot':
                    self.changes['heatmap_action'] = flatten[at]
                elif k[1] == 'General':
                    if len(k) > 5 and k[5] in ['dpi-save', 'save-format'] or k[3] == 'toolbar':
                        self.changes['figure'] = flatten[at]
                    else:
                        self.changes['bar_action'] = flatten[at]
                        self.changes['line_ie_action'] = flatten[at]
                        self.changes['line_action'] = flatten[at]
                        self.changes['heatmap_action'] = flatten[at]
                else:
                    self.changes['visualization_action'] = flatten[at]
        self.update(osett)

    def get_settings(self):
        flatten = flatten_dict(self)

        items = [
            (tuple(ik for ik in k if ik not in ['children', 'value']), v)
            for k, v in flatten.items()
            if 'action_type' not in k
        ]
        return dict(items)

    def set_as_default(self):
        value2deafult(self)


class CorrChartSettings(dict):
    D = 1  # Drawable
    U = 2  # Updatable
    R = 3  # re-plot
    default = {
        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'parameters', 'value': None, 'default': None,
        'children': {
            'General': {
                'type': 'group', 'enabled': True, 'expanded': True, 'name': 'General', 'value': None, 'default':
                    None,
                'children': {
                    'theme': {
                        'type': 'list', 'enabled': True, 'expanded': True, 'name': 'theme', 'value': 'darkgrid',
                        'values': ['darkgrid', '---0', 'whitegrid', 'dark', 'white', 'ticks'], 'default': 'darkgrid',
                        'action_type': R, 'tip': 'changes visual theme for all plots'},
                    'toolbar': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'toolbar', 'value': False,
                                'default': False, 'action_type': U, 'tip': 'shows toolbar in all graphs'},
                    'figure-format': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'figure-format', 'value': None,
                        'default': None,
                        'children': {
                            'dpi-plot': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'dpi-plot',
                                         'value': 100, 'step': 10, 'accelerated': True, 'limits': (50, 300),
                                         'default': 100, 'action_type': R, 'tip': 'set dpi for all plots'},
                            'dpi-save': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'dpi-save',
                                         'value': 300, 'step': 10, 'accelerated': True, 'limits': (50, 1200),
                                         'default': 300, 'action_type': U, 'tip': 'dpi used for saved figures'},
                            'save-format': {
                                'type': 'list', 'enabled': True, 'expanded': True, 'name': 'save-format',
                                'value': 'svg', 'values': ['svg', 'png', 'tiff', 'pdf', 'eps', 'jpg', 'jpeg', 'pgf',
                                                           'ps', 'raw', 'rgba', 'svgz', 'tif'],
                                'default': 'svg', 'action_type': U},
                        }}}},
            'Regression Plot': {
                'type': 'group', 'enabled': True, 'expanded': False, 'name': 'Regression Plot', 'value': None,
                'default': None,
                'children': {
                    'line-width': {'type': 'float', 'enabled': True, 'expanded': True, 'name': 'line-width',
                                   'value': 1.0, 'step': 0.1, 'limits': (0.1, 1.5), 'accelerated': True,
                                   'default': 1.0, 'action_type': R},
                    'line-color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'line-color',
                                   'value': [0, 0, 0, 255], 'default': [0, 0, 0, 255], 'action_type': R},
                    'Conf. Interval': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'Conf. Interval',
                                       'value': 95, 'step': 1, 'limits': (0, 100), 'accelerated': True,
                                       'default': 95, 'action_type': R},
                    'pearson': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'pearson', 'value': True,
                                'default': True, 'action_type': R},
                    'spearman': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'spearman', 'value': True,
                                 'default': True, 'action_type': R},
                    'equation': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'equation', 'value': True,
                                 'default': True, 'action_type': R},
                    'Scatter': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'Scatter', 'value': None,
                        'default': None,
                        'children': {
                            'color': {'type': 'color', 'enabled': True, 'expanded': True, 'name': 'color',
                                      'value': [44, 105, 176, 255], 'default': [44, 105, 176, 255], 'action_type': R},
                            'alpha': {'type': 'float', 'enabled': True, 'expanded': True, 'name': 'alpha',
                                      'value': 0.15, 'step': 0.01, 'default': 0.15, 'limits': (0, 1), 'action_type': R},
                            'marker-size': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'marker-size',
                                            'value': 10, 'step': 1, 'default': 10, 'limits': (2, 100),
                                            'accelerated': True, 'action_type': R},
                            'error-line': {'type': 'group', 'enabled': True, 'expanded': True, 'name': 'error-line',
                                           'value': None, 'default': None,
                                           'children': {
                                               'show': {'type': 'bool', 'enabled': True, 'expanded': True,
                                                        'name': 'show', 'value': False, 'default': False,
                                                        'action_type': R},
                                               'width': {'type': 'float', 'enabled': True, 'expanded': True,
                                                         'name': 'width',
                                                         'value': 0.7, 'step': 0.1, 'limits': (0.1, 1.5),
                                                         'accelerated': True, 'default': 0.7, 'action_type': R},
                                               'color': {'type': 'color', 'enabled': True, 'expanded': True,
                                                         'name': 'color', 'value': [0, 0, 0, 255],
                                                         'default': [0, 0, 0, 255], 'action_type': R},
                                               'cap-size': {'type': 'int', 'enabled': True, 'expanded': True,
                                                            'name': 'cap-size', 'value': 0, 'step': 1,
                                                            'limits': (0, 50),
                                                            'accelerated': True, 'default': 0, 'action_type': R},
                                               'representation': {'type': 'list', 'enabled': True, 'expanded': True,
                                                                  'name': 'representation', 'value': 'SD',
                                                                  'values': ['SD', 'SEM'], 'default': 'SD',
                                                                  'action_type': R,
                                                                  'tip': "The metric representing the error line. SD, the "
                                                                         "Standard Deviation, describes the data distribution "
                                                                         "around the mean, while SEM, the Standard Error of "
                                                                         "the Mean, describes how representative the mean is "
                                                                         "of the population."}}}}},
                    'Distribution': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'Distribution', 'value': None,
                        'default': None,
                        'children': {
                            'show': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'show', 'value': True,
                                     'default': True, 'action_type': R},
                            'type': {'type': 'list', 'enabled': True, 'expanded': True, 'name': 'type',
                                      'value': 'hist', 'values': ['hist', 'kde', 'violin', 'box'],
                                      'default': 'hist', 'action_type': R},
                            'Histogram': {
                                'type': 'group', 'enabled': True, 'expanded': True, 'name': 'Histogram',
                                'value': None, 'default': None,
                                'children': {
                                    'fill-bars': {'type': 'bool', 'enabled': True, 'expanded': True,
                                                  'name': 'fill-bars', 'value': False, 'default': False,
                                                  'action_type': R},
                                    'kde': {'type': 'bool', 'enabled': True, 'expanded': True, 'name': 'kde',
                                            'value': True, 'default': True, 'action_type': R},
                                    'element': {'type': 'list', 'enabled': True, 'expanded': True,
                                                'name': 'element', 'value': 'bars', 'values': ['step', 'bars', 'poly'],
                                                'default': 'bars', 'action_type': R}}}}},
                    'fontsize': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'fontsize', 'value': None,
                        'default': None,
                        'children': {
                            'x-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'y-ticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-ticks',
                                        'value': 10, 'default': 10, 'action_type': D},
                            'x-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'x-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'y-label': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'y-label',
                                        'value': 12, 'default': 12, 'action_type': D},
                            'title': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'title',
                                      'value': 14, 'default': 14, 'action_type': D},
                            'suptitle': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'suptitle',
                                         'value': 12, 'default': 12, 'action_type': D},
                            'legend': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'legend',
                                       'value': 9, 'default': 9, 'action_type': D}}},
                    'axes': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'axes', 'value': None,
                        'default': None,
                        'children': {
                            'num-xticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'num-xticks',
                                           'value': 10, 'default': 10, 'action_type': R},
                            'num-yticks': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'num-yticks',
                                           'value': 10, 'default': 10, 'action_type': R}}},
                    'figure': {
                        'type': 'group', 'enabled': True, 'expanded': True, 'name': 'figure', 'value': None,
                        'default': None,
                        'children': {
                            'width': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'width', 'value': 6,
                                      'default': 6, 'action_type': R},
                            'height': {'type': 'int', 'enabled': True, 'expanded': True, 'name': 'height',
                                       'value': 6, 'default': 6, 'action_type': R}}}}
            },
        }}

    def __init__(self):
        super(CorrChartSettings, self).__init__()

        self.update(self.default)

        self.changes = 0

    def return_default(self):
        self.update(self.default)

    def is_changed(self, osett) -> bool:
        flatten = flatten_dict(self)
        f_osett = flatten_dict(osett)
        return flatten != f_osett

    def get_changes(self, osett) -> None:

        flatten = flatten_dict(self)
        f_osett = flatten_dict(osett)
        if flatten == f_osett:
            return
        for k, v in flatten.items():
            if k[-1] == 'value' and v != f_osett[k]:
                at = tuple(list(k[:-1]) + ['action_type'])
                self.changes = flatten[at]
        self.update(osett)

    def get_settings(self):
        flatten = flatten_dict(self)

        items = [
            (tuple(ik for ik in k if ik not in ['children', 'value']), v)
            for k, v in flatten.items()
            if 'action_type' not in k
        ]
        return dict(items)


class Palette(LinearSegmentedColormap):
    def __init__(self, name, clist, ptype):
        super(Palette, self).__init__(name, clist)
        self.name = name
        self.colors = [list(map(lambda x: x / 255, color)) for color in clist]
        self.ptype = ptype
        self.colormap = self.from_list(self.name, self.colors)


color_store = {
    'cartoon': {
        'sequential': {
        },
        'diverging': dict(
            ArmyRose_5=[
                [121, 130, 52],
                [208, 211, 162],
                [253, 251, 228],
                [240, 198, 195],
                [212, 103, 128]],
            ArmyRose_7=[
                [121, 130, 52],
                [163, 173, 98],
                [208, 211, 162],
                [253, 251, 228],
                [240, 198, 195],
                [223, 145, 163],
                [212, 103, 128]],
            Geyser_5=[
                [0, 128, 128],
                [180, 200, 168],
                [246, 237, 189],
                [237, 187, 138],
                [202, 86, 44]],
            Geyser_7=[
                [0, 128, 128],
                [112, 164, 148],
                [180, 200, 168],
                [246, 237, 189],
                [237, 187, 138],
                [222, 138, 90],
                [202, 86, 44]],
            Temps_5=[
                [0, 147, 146],
                [156, 203, 134],
                [233, 226, 156],
                [238, 180, 121],
                [207, 89, 126]],
            Temps_7=[
                [0, 147, 146],
                [57, 177, 133],
                [156, 203, 134],
                [233, 226, 156],
                [238, 180, 121],
                [232, 132, 113],
                [207, 89, 126]],
            TealRose_5=[
                [0, 147, 146],
                [177, 199, 179],
                [241, 234, 200],
                [229, 185, 173],
                [208, 88, 126]],
            TealRose_7=[
                [0, 147, 146],
                [114, 170, 161],
                [177, 199, 179],
                [241, 234, 200],
                [229, 185, 173],
                [217, 137, 148],
                [208, 88, 126]],
            Tropic_5=[
                [0, 155, 158],
                [167, 211, 212],
                [241, 241, 241],
                [228, 193, 217],
                [199, 93, 171]],
            Tropic_7=[
                [0, 155, 158],
                [66, 183, 185],
                [167, 211, 212],
                [241, 241, 241],
                [228, 193, 217],
                [214, 145, 193],
                [199, 93, 171]],
        ),
        'Qualitative': {
        }},
    'cmocean': {
        'diverging': dict(
            Balance_5=[
                [24, 28, 67],
                [56, 136, 186],
                [241, 236, 235],
                [192, 90, 60],
                [60, 9, 18]],
            Balance_7=[
                [24, 28, 67],
                [12, 94, 190],
                [117, 170, 190],
                [241, 236, 235],
                [208, 139, 115],
                [167, 36, 36],
                [60, 9, 18]],
            Curl_5=[
                [21, 29, 68],
                [44, 148, 127],
                [254, 246, 245],
                [196, 90, 97],
                [52, 13, 53]],
            Curl_7=[
                [21, 29, 68],
                [21, 109, 115],
                [125, 179, 144],
                [254, 246, 245],
                [219, 140, 119],
                [157, 48, 96],
                [52, 13, 53]],
            Delta_5=[
                [17, 32, 64],
                [51, 145, 169],
                [255, 253, 205],
                [97, 146, 11],
                [23, 35, 19]],
            Delta_7=[
                [17, 32, 64],
                [28, 103, 160],
                [108, 181, 179],
                [255, 253, 205],
                [170, 172, 32],
                [24, 115, 40],
                [23, 35, 19]])
    },
    'colorbrewer': {
        'diverging': dict(
            BrBG_5=[[166, 97, 26],
                    [223, 194, 125],
                    [245, 245, 245],
                    [128, 205, 193],
                    [1, 133, 113]],
            BrBG_7=[[140, 81, 10],
                    [216, 179, 101],
                    [246, 232, 195],
                    [245, 245, 245],
                    [199, 234, 229],
                    [90, 180, 172],
                    [1, 102, 94]],
            PRGn_5=[[123, 50, 148],
                    [194, 165, 207],
                    [247, 247, 247],
                    [166, 219, 160],
                    [0, 136, 55]],
            PRGn_7=[[118, 42, 131],
                    [175, 141, 195],
                    [231, 212, 232],
                    [247, 247, 247],
                    [217, 240, 211],
                    [127, 191, 123],
                    [27, 120, 55]],
            PiYG_5=[[208, 28, 139],
                    [241, 182, 218],
                    [247, 247, 247],
                    [184, 225, 134],
                    [77, 172, 38]],
            PiYG_7=[[197, 27, 125],
                    [233, 163, 201],
                    [253, 224, 239],
                    [247, 247, 247],
                    [230, 245, 208],
                    [161, 215, 106],
                    [77, 146, 33]],
            PuOr_5=[[230, 97, 1],
                    [253, 184, 99],
                    [247, 247, 247],
                    [178, 171, 210],
                    [94, 60, 153]],
            PuOr_7=[[179, 88, 6],
                    [241, 163, 64],
                    [254, 224, 182],
                    [247, 247, 247],
                    [216, 218, 235],
                    [153, 142, 195],
                    [84, 39, 136]],
            RdBu_5=[[202, 0, 32],
                    [244, 165, 130],
                    [247, 247, 247],
                    [146, 197, 222],
                    [5, 113, 176]],
            RdBu_7=[[178, 24, 43],
                    [239, 138, 98],
                    [253, 219, 199],
                    [247, 247, 247],
                    [209, 229, 240],
                    [103, 169, 207],
                    [33, 102, 172]],
            RdGy_5=[[202, 0, 32],
                    [244, 165, 130],
                    [255, 255, 255],
                    [186, 186, 186],
                    [64, 64, 64]],
            RdGy_7=[[178, 24, 43],
                    [239, 138, 98],
                    [253, 219, 199],
                    [255, 255, 255],
                    [224, 224, 224],
                    [153, 153, 153],
                    [77, 77, 77]])
    }
}

seismic_colors = [list(map(lambda x: x * 255, c)) for c in datad['seismic']]


class Palettes:
    # -- BEGIN -- Diverging palettes
    # Cartoon
    ArmyRose_5 = Palette('Armyrose_5', color_store['cartoon']['diverging']['ArmyRose_5'], 'diverging')
    ArmyRose_7 = Palette('Armyrose_7', color_store['cartoon']['diverging']['ArmyRose_7'], 'diverging')
    Geyser_5 = Palette('Geyser_5', color_store['cartoon']['diverging']['Geyser_5'], 'diverging')
    Geyser_7 = Palette('Geyser_7', color_store['cartoon']['diverging']['Geyser_7'], 'diverging')
    Temps_5 = Palette('Temps_5', color_store['cartoon']['diverging']['Temps_5'], 'diverging')
    Temps_7 = Palette('Temps_7', color_store['cartoon']['diverging']['Temps_7'], 'diverging')
    TealRose_5 = Palette('TealRose_5', color_store['cartoon']['diverging']['TealRose_5'], 'diverging')
    TealRose_7 = Palette('TealRose_7', color_store['cartoon']['diverging']['TealRose_7'], 'diverging')
    Tropic_5 = Palette('Tropic_5', color_store['cartoon']['diverging']['Tropic_5'], 'diverging')
    Tropic_7 = Palette('Tropic_7', color_store['cartoon']['diverging']['Tropic_7'], 'diverging')
    # cmocean
    Balance_5 = Palette('Balance_5', color_store['cmocean']['diverging']['Balance_5'], 'diverging')
    Balance_7 = Palette('Balance_7', color_store['cmocean']['diverging']['Balance_7'], 'diverging')
    Curl_5 = Palette('Curl_5', color_store['cmocean']['diverging']['Curl_5'], 'diverging')
    Curl_7 = Palette('Curl_7', color_store['cmocean']['diverging']['Curl_7'], 'diverging')
    Delta_5 = Palette('Delta_5', color_store['cmocean']['diverging']['Delta_5'], 'diverging')
    Delta_7 = Palette('Delta_7', color_store['cmocean']['diverging']['Delta_7'], 'diverging')
    # colorbrew
    BrBG_5 = Palette('BrBG_5', color_store['colorbrewer']['diverging']['BrBG_5'], 'diverging')
    BrBG_7 = Palette('BrBG_7', color_store['colorbrewer']['diverging']['BrBG_7'], 'diverging')
    PRGn_5 = Palette('PRGn_5', color_store['colorbrewer']['diverging']['PRGn_5'], 'diverging')
    PRGn_7 = Palette('PRGn_7', color_store['colorbrewer']['diverging']['PRGn_7'], 'diverging')
    PiYG_5 = Palette('PiYG_5', color_store['colorbrewer']['diverging']['PiYG_5'], 'diverging')
    PiYG_7 = Palette('PiYG_7', color_store['colorbrewer']['diverging']['PiYG_7'], 'diverging')
    PuOr_5 = Palette('PuOr_5', color_store['colorbrewer']['diverging']['PuOr_5'], 'diverging')
    PuOr_7 = Palette('PuOr_7', color_store['colorbrewer']['diverging']['PuOr_7'], 'diverging')
    RdBu_5 = Palette('RdBu_5', color_store['colorbrewer']['diverging']['RdBu_5'], 'diverging')
    RdBu_7 = Palette('RdBu_7', color_store['colorbrewer']['diverging']['RdBu_7'], 'diverging')
    RdGy_5 = Palette('RdGy_5', color_store['colorbrewer']['diverging']['RdGy_5'], 'diverging')
    RdGy_7 = Palette('RdGy_7', color_store['colorbrewer']['diverging']['RdGy_7'], 'diverging')
    # matplolib
    seismic = Palette('seismic', seismic_colors, 'diverging')
    RdYlBu = 'RdYlBu'
    RdYlGn = 'RdYlGn'
    Spectral = 'Spectral'
    coolwarm = 'coolwarm'
    # -- END -- Diverging palettes

    # -- BEGIN -- Qualitative colormaps
    # seaborn
    deep = 'deep'
    muted = 'muted'
    bright = 'bright'
    pastel = 'pastel'
    dark = 'dark'
    colorblind = 'colorblind'
    # matplotlib
    Pastel1 = 'Pastel1'
    Pastel2 = 'Pastel2'
    Paired = 'Paired'
    Accent = 'Accent'
    Dark2 = 'Dark2'
    Set1 = 'Set1'
    Set2 = 'Set2'
    Set3 = 'Set3'
    tab10 = 'tab10'
    tab20 = 'tab20'
    tab20b = 'tab20b'
    tab20c = 'tab20c'
    husl = 'husl'
    hls = 'hls'

    # -- END -- Qualitative colormaps

    @classmethod
    def get_colormap(cls, name):
        palette = getattr(cls, name)
        if isinstance(palette, Palette):
            return palette.colormap
        else:
            return palette

    @classmethod
    def get_palette(cls, name):
        return getattr(cls, name)
