# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA                  #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

from GMXMMPBSA.gui import GMX_MMPBSA_GUI
from GMXMMPBSA.commandlineparser import guiparser
from PyQt5.QtWidgets import QApplication
import sys
from pathlib import Path

if __name__ == '__main__':
    app = QApplication(sys.argv)
    parser = guiparser.parse_args(sys.argv[1:])
    path = Path(parser.path).absolute()
    if not path.exists():
        print('Path not found')
        sys.exit(1)
    app.setApplicationName('GMX-MMPBSA-GUI')
    w = GMX_MMPBSA_GUI(path.as_posix())
    w.show()
    sys.exit(app.exec())
