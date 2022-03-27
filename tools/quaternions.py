import numpy as np


class Quaternion:
    """
    Quaternion class
    """

    def __init__(self, q0, q1, q2, q3):
        self.q0 = q0
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3

    @classmethod
    def from_angles(cls, theta, phi, psi):
        """
        Create quaternion from rotation angles:
        theta: angle of rotation
        phi: azimuthal angle of axis of rotation
        psi: polar angle of axis of rotation
        """
        return cls(
            np.cos(theta / 2),
            np.sin(theta / 2) * np.cos(phi) * np.sin(psi),
            np.sin(theta / 2) * np.sin(phi) * np.sin(psi),
            np.sin(theta / 2) * np.cos(psi),
        )

    @classmethod
    def from_array(cls, a):
        return cls(a[0], a[1], a[2], a[3])

    @classmethod
    def from_scalar_vector(cls, q0, q):
        return cls(q0, q[0], q[1], q[2])

    @property
    def q(self):
        return np.array([self.q1, self.q2, self.q3])

    @q.setter
    def q(self, new_q):
        self.q1 = new_q[0]
        self.q2 = new_q[1]
        self.q3 = new_q[2]

    @property
    def abs(self):
        return q0 ** 2 + np.dot(self.q, self.q)

    def __str__(self):
        return f"({self.q0}, {self.q1}, {self.q2}, {self.q3})"

    def __repr__(self):
        return f"quaternions.Quaternion({self.q0}, {self.q1}, {self.q2}, {self.q3})"

    def __getitem__(self, arg):
        return self.to_array()[arg]

    def __add__(self, other):
        return Quaternion.from_array(self.to_array() + other.to_array())

    def __mul__(self, other):
        if type(other) == Quaternion:
            return Quaternion.from_scalar_vector(
                self.q0 * other.q0 - np.dot(self.q, other.q),
                self.q0 * other.q + self.q * other.q0 + np.cross(self.q, other.q),
            )
        else:
            return Quaternion.from_array(other * self.to_array())

    def to_array(self):
        return np.array([self.q0, self.q1, self.q2, self.q3])

    def to_rotation_matrix(self):
        (q0, q1, q2, q3) = tuple(self.to_array())
        return np.array(
            [
                [2 * (q0 ** 2 + q1 ** 2) - 1, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2),],
                [2 * (q1 * q2 + q0 * q3), 2 * (q0 ** 2 + q2 ** 2) - 1, 2 * (q2 * q3 - q0 * q1),],
                [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 2 * (q0 ** 2 + q3 ** 2) - 1,],
            ]
        )
