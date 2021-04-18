import numpy as np


class EoMSolver:
    """
    Equation of motion solver
    """
    def __init__(self):
        self.m = 0
        self.k = 0
        self.Z = 0
        self.u0 = 0
        self.v0 = 0
        self.f = 0


    def equationOfMotion(self, v, u, t):
        return (self.f(t) - self.c * v - self.k * u) / self.m


    def mass(self, m):
        if isinstance(m, (int, float)):
            self.m = m


    def stiffness(self, k):
        if isinstance(k, (int, float)):
            self.k = k


    def damping(self, Z):
        if isinstance(Z, (int, float)):
            self.Z = Z


    def force(self, f):
        self.f = f


    def initialCondition(self, u0, v0):
        if isinstance(u0, (int, float)):
            self.u0 = u0
        if isinstance(v0, (int, float)):
            self.v0 = v0


    def solve(self, time_max, step_size):
        # Damping
        self.c = self.Z / 100 * 2 * (self.k * self.m) ** 0.5

        # time steps
        time_steps = int(time_max / step_size) + 1
        self.t = np.linspace(0, time_max, time_steps)

        # size
        size = self.t.size

        # zeroing
        self.u = np.zeros([size])
        self.v = np.zeros([size])

        # initial conditions
        self.u[0] = self.u0
        self.v[0] = self.v0

        # Solving
        for self.i in range(size - 1):
            # print(self.i + 1)
            self.v[self.i + 1], self.u[self.i + 1] = self.nextStep()

        return self.v, self.u, self.t

    def nextStep(self):
        pass


class RungeKutta4th(EoMSolver):
    def nextStep(self):
        v, u, f, i, t = self.v, self.u, self.f, self.i, self.t
        dt = t[i + 1] - t[i]

        K1u = dt * v[i]
        K1v = dt * self.equationOfMotion(v[i], u[i], t[i])

        K2u = dt * (v[i] + K1v / 2)
        K2v = dt * self.equationOfMotion(v[i] + K1v / 2, u[i] + K1u / 2, t[i] + dt / 2)

        K3u = dt * (v[i] + K2v / 2)
        K3v = dt * self.equationOfMotion(v[i] + K2v / 2, u[i] + K2u / 2, t[i] + dt / 2)

        K4u = dt * (v[i] + K3v)
        K4v = dt * self.equationOfMotion(v[i] + K3v, u[i] + K3u, t[i] + dt)

        return v[i] + (K1v + 2 * K2v + 2 * K3v + K4v) / 6, u[i] + (K1u + 2 * K2u + 2 * K3u + K4u) / 6


class ForwardEuler(EoMSolver):
    def nextStep(self):
        v, u, f, i, t = self.v, self.u, self.f, self.i, self.t
        dt = t[i + 1] - t[i]
        return v[i] + dt * self.equationOfMotion(v[i], u[i], t[i]), u[i] + dt * v[i]


class LeapFrog(EoMSolver):
    def nextStep(self):
        v, u, f, i, t = self.v, self.u, self.f, self.i, self.t
        dt = t[i + 1] - t[i]

        # Starting velocity
        if i == 0:
            self.v_lastHalfStep = v[i] - self.equationOfMotion(v[i], u[i], t[i]) * dt / 2

        v_nextHalfStep = self.v_lastHalfStep + self.equationOfMotion(v[i], u[i], t[i]) * dt
        u_next = u[i] + v_nextHalfStep * dt
        v_next = (v_nextHalfStep + self.v_lastHalfStep) / 2
        self.v_lastHalfStep = v_nextHalfStep
        return v_next, u_next


class Newmark(EoMSolver):
    def nextStep(self):
        v, u, f, i, t, k, m, c = self.v, self.u, self.f, self.i, self.t, self.k, self.m, self.c
        dt = t[i + 1] - t[i]
        gamma = 0.5
        beta = 1.03
        k_prime = k + gamma * c / beta / dt + m / beta / (dt**2)
        f_prime = f(t[i + 1]) + (m / beta / (dt**2) + gamma * c / beta / dt) * u[i] + (m / beta / dt + (gamma / beta - 1) * c) * v[i] + ((0.5 * beta - 1) * m + dt * (gamma / 2 / beta - 1) * c) * self.equationOfMotion(v[i], u[i], t[i])
        u_next = f_prime / k_prime
        return gamma / beta / dt * (u_next - u[i]) + (1 - gamma / beta) * v[i] + dt * (1 - gamma / 2 / beta) * self.equationOfMotion(v[i], u[i], t[i]), u_next


if __name__ == '__main__':
    from matplotlib import pyplot as plt


    def f(t):
        f = 100*np.sin(1.2 * t) if t < 90 else 0
        return f


    for solver in [RungeKutta4th, Newmark, ForwardEuler, LeapFrog]:
        case = solver()
        case.stiffness(100)
        case.mass(100)
        case.damping(1)
        case.force(f)
        case.initialCondition(1, 0)
        v, u, t = case.solve(50, 0.01)
        plt.plot(t, u, label=solver.__name__)

    plt.legend()
    plt.show()