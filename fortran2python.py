
from cc3d import CompuCellSetup
        

from fortran2pythonSteppables import fortran2pythonSteppable

CompuCellSetup.register_steppable(steppable=fortran2pythonSteppable(frequency=1))


CompuCellSetup.run()
