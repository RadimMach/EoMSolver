import sys
from EoMSolver import *
from PyQt5.QtWidgets import *
from PyQt5.uic import loadUi
from PyQt5.QtGui import QDoubleValidator, QValidator, QRegExpValidator
from PyQt5.QtCore import QRegExp, QThread, pyqtSignal, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure


class MyWindow(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi()

    def setupUi(self):
        # Load GUI
        self.ui = loadUi("Dialog.ui", self)

        # Create Matplotlib plot
        self.canvas = MplCanvas(self)
        self.ui.sdof_dialog.addWidget(self.canvas, 1, 1, 1, 1)
        # self.addToolBar(NavigationToolbar(self.canvas, self))

        # Generate central widget
        self.setCentralWidget(self.centralWidget)

        # Pushbutton
        self.pushButton.clicked.connect(self.checkInput)

        # Progress bar
        self.progressBar.hide()

        # Validate input
        validator = QDoubleValidator(0.001, 1000000, 3)
        validator_2 = QDoubleValidator(0.000, 1000000, 3)
        validator_3 = QDoubleValidator(-1000000, 1000000, 3)

        self.validating(self.ui.mass_lineEdit, validator)
        self.validating(self.ui.stiffness_lineEdit, validator)
        self.validating(self.ui.damping_lineEdit, validator_2)
        self.validating(self.ui.initial_displacement_lineEdit, validator_3)
        self.validating(self.ui.initial_velocity_lineEdit, validator_3)
        self.validating(self.ui.end_time_lineEdit, validator)
        self.validating(self.ui.time_step_lineEdit, validator)

        # Connect checkboxes
        self.ui.displacement_radioButton.toggled.connect(self.toPlot)
        self.ui.velocity_radioButton.toggled.connect(self.toPlot)
        self.ui.RungeKutta4th_checkBox.stateChanged.connect(self.toPlot)
        self.ui.ForwardEuler_checkBox.stateChanged.connect(self.toPlot)
        self.ui.LeapFrog_checkBox.stateChanged.connect(self.toPlot)
        self.ui.Newmark_checkBox.stateChanged.connect(self.toPlot)


        self.ui.mass_lineEdit.setText("100")
        self.ui.stiffness_lineEdit.setText("100")
        self.ui.initial_displacement_lineEdit.setText("1")
        self.ui.end_time_lineEdit.setText("100")
        self.ui.time_step_lineEdit.setText("0.01")
        self.ui.excitation_force_plainTextEdit.setPlainText("100*sin(1.2*t) if t < 60 else 0")

        # Zeroing solution
        self.solutions = None


    def errorMessage(self, message):
        error_dialog = QMessageBox()
        error_dialog.setIcon(QMessageBox.Warning)
        error_dialog.setWindowTitle("Warning!")
        error_dialog.setText(message)
        error_dialog.exec_()

    def validating(self, toValidate, validationRule):

        def checkIfValid():
            if validationRule.validate(toValidate.text(), 1)[0] == QValidator.Acceptable:
                pass
            # to ignore second signal - bug
            elif not toValidate.isModified():
                toValidate.setModified(True)
            # if non-valid input bring Error message
            else:
                toValidate.setModified(False)
                toValidate.setText("")
                toValidate.setFocus()

        # Validate object input
        toValidate.editingFinished.connect(checkIfValid)

    def checkInput(self):
        # Check input
        validLineEdit = True
        validPLainText = True
        stringCallable = True
        methodToSolve = True

        # Get children of verticalLayout
        for count in range(self.ui.parameters_verticalLayout.count()):
            child = self.ui.parameters_verticalLayout.itemAt(count).widget()

            # Check if specified lineEdit are differ then zero
            if child.objectName() in ["mass_lineEdit", "stiffness_lineEdit", "end_time_lineEdit", "time_step_lineEdit"]:
                if child.text() == "":
                    validLineEdit = False

            # Check if excitation force is valid
            elif child.objectName() == "excitation_force_plainTextEdit":
                excitationForce = child.toPlainText()
                validator_function = QRegExpValidator(QRegExp(r"[\d|t|sin|cos| |+|-|*|/|.|(|)|=|<|>|if|else|elif|:]*"))

                # If string contains valid characters
                if validator_function.validate(excitationForce, 2)[0] == QValidator.Acceptable:
                    excitationForce = excitationForce.replace("cos", "np.cos")
                    excitationForce = excitationForce.replace("sin", "np.sin")
                    # Convert string to callable function
                    self.f = lambda t: eval(excitationForce)

                    # Test if callable string is valid
                    try:
                        self.f(0)
                    except:
                        stringCallable = False

                # If string contains also non-valid characters
                else:
                    validPLainText = False
        # Any method to solve
        methodChecked = [self.ui.solve_ForwardEuler_checkBox.isChecked(),
                         self.ui.solve_RungeKutta4th_checkBox.isChecked(),
                         self.ui.solve_LeapFrog_checkBox.isChecked(),
                         self.ui.solve_Newmark_checkBox.isChecked()]
        if not any(methodChecked):
            methodToSolve = False


        # If error show message else solve
        if not validLineEdit:
            self.errorMessage("Mass, Stiffness, End time and Time step must be bigger than 0")
        elif float(self.ui.end_time_lineEdit.text()) <= float(self.ui.time_step_lineEdit.text()):
            self.errorMessage("Time step must be smaller than End time")
        elif not validPLainText:
            self.errorMessage("Excitation force contains forbidden characters")
        elif not stringCallable:
            self.errorMessage("Check excitation force, it is not written correctly")
        elif not methodToSolve:
            self.errorMessage("Choose at least one method to solve")
        else:
            self.solveInput()


    def solveInput(self):
        # Get inputs
        m = float(self.ui.mass_lineEdit.text())
        k = float(self.ui.stiffness_lineEdit.text())
        c = float(self.ui.damping_lineEdit.text())
        u = float(self.ui.initial_displacement_lineEdit.text())
        v = float(self.ui.initial_velocity_lineEdit.text())
        t = float(self.ui.end_time_lineEdit.text())
        dt = float(self.ui.time_step_lineEdit.text())

        # Get methods to solve
        methods = []

        if self.ui.solve_ForwardEuler_checkBox.isChecked():
            methods.append(ForwardEuler)
        if self.ui.solve_RungeKutta4th_checkBox.isChecked():
            methods.append(RungeKutta4th)
        if self.ui.solve_LeapFrog_checkBox.isChecked():
            methods.append(LeapFrog)
        if self.ui.solve_Newmark_checkBox.isChecked():
            methods.append(Newmark)

        # Disable input
        self.inputEnabled(False)

        # Show Progress bar
        self.ui.progressBar.show()
        self.ui.progressBar.setValue(0)

        # Delete solutions
        self.solutions = None

        # Start thread to solve input
        self.solver = SolverThread(m, k, c, u, v, t, dt, self.f, methods)
        self.solver.start()
        self.solver.updateProgress.connect(self.updateProgressBar)
        self.solver.solverFinished.connect(self.finishSolver)
        self.solver.finished.connect(self.toPlot)


    def toPlot(self):

        # If there is any solution
        if self.solutions is not None:
            # Clear figure
            self.canvas.clearFigure()

            # Methods checkboxes
            methodsChecked = {"RungeKutta4th": self.ui.RungeKutta4th_checkBox.isChecked(),
                              "ForwardEuler":  self.ui.ForwardEuler_checkBox.isChecked(),
                              "LeapFrog":      self.ui.LeapFrog_checkBox.isChecked(),
                              "Newmark":       self.ui.Newmark_checkBox.isChecked()}

            quantityChecked = {"Displacement": self.ui.displacement_radioButton.isChecked(),
                               "Velocity":     self.ui.velocity_radioButton.isChecked()}

            legend = []
            title = ""

            for (method, quantity), solution in self.solutions.items():
                # If method checkbox and data checkbox are checked
                if methodsChecked[method] and quantityChecked[quantity]:
                    # Plot data
                    self.canvas.setPlot(self.t, solution)
                    # Legend and Title
                    legend.append(method)
                    title = quantity

            # draw plots with settings
            self.canvas.drawPlots(legend, title)


    def updateProgressBar(self, progress):
        self.ui.progressBar.setValue(progress)


    def inputEnabled(self, state=True):
        self.pushButton.setEnabled(state)
        self.ui.solve_ForwardEuler_checkBox.setEnabled(state)
        self.ui.solve_LeapFrog_checkBox.setEnabled(state)
        self.ui.solve_Newmark_checkBox.setEnabled(state)
        self.ui.solve_RungeKutta4th_checkBox.setEnabled(state)
        self.ui.excitation_force_plainTextEdit.setEnabled(state)
        self.ui.time_step_lineEdit.setEnabled(state)
        self.ui.end_time_lineEdit.setEnabled(state)
        self.ui.initial_velocity_lineEdit.setEnabled(state)
        self.ui.initial_displacement_lineEdit.setEnabled(state)
        self.ui.damping_lineEdit.setEnabled(state)
        self.ui.stiffness_lineEdit.setEnabled(state)
        self.ui.mass_lineEdit.setEnabled(state)


    def finishSolver(self, solverOutput):
        # Get solutions
        self.solutions, self.t = solverOutput

        # Enable input
        self.inputEnabled()

        #Hide progress bar
        self.ui.progressBar.hide()


class SolverThread(QThread):
    """
    Class for solving in another thread
    """
    updateProgress = pyqtSignal(int)
    solverFinished = pyqtSignal(list)

    def __init__(self, m, k, c, u, v, t, dt, f, methods, parent=None):
        QThread.__init__(self, parent)
        self.m = m
        self.k = k
        self.c = c
        self.u = u
        self.v = v
        self.t = t
        self.dt = dt
        self.f = f
        self.methods = methods


    def run(self):
        # Zeroing matrixes
        solutions = {}

        for count, method in enumerate(self.methods):
            # Set input to solver
            solver = method()
            solver.mass(self.m)
            solver.stiffness(self.k)
            solver.damping(self.c)
            solver.initialCondition(self.u, self.v)
            solver.force(self.f)

            # Solve
            vMethod, uMethod, time = solver.solve(self.t, self.dt)

            # Solution to dictionary
            solutions[(method.__name__, "Velocity")] = vMethod
            solutions[(method.__name__, "Displacement")] = uMethod

            #Update progress
            self.updateProgress.emit(int((count+1)/len(self.methods)*100))

        # solution to list
        self.solverFinished.emit([solutions, time])


class MplCanvas(FigureCanvas):
    """
    Class for integration MatPlotLib to QWidget
    """
    def __init__(self, parent=None):
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        super(MplCanvas, self).__init__(self.figure)
        self.figure.tight_layout()


    def clearFigure(self):
        self.axes.clear()


    def setPlot(self, t, variable):
        self.axes.plot(t, variable)


    def drawPlots(self, legend, title):
        self.axes.legend((legend), loc='upper right')
        self.axes.set_title(title)
        self.figure.tight_layout()
        self.draw()

app = QApplication([])
window = MyWindow()
window.show()
sys.exit(app.exec())
