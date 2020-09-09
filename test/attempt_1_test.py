import unittest
import numpy as np
from src.attempt_1 import acceleration#, runkut

M0 = 2 * 10 ** 30  ##solar mass
m_e = 6 * 10 ** 24  ## earth mass
R = 1.496 * 10 ** 8  ## astronomical unit
G = 6.67 * 10 ** (-11)
V = 3 * 10 ** 4  ## earth avg orbital speed


class TestMethods(unittest.TestCase):

    def test_acceleration(self):
        dummy_r = np.array([3, 4, 12])
        dummy_m2 = 10
        #norm = 13
        #numerator = - 6.67 * 10 ** (-10)
        #norm_sq = 169
        #multiplier = -3.946745562130177 * 10 ** -12
        expected_dvdt = np.array([-9.10787437e-13, -1.21438325e-12, -3.64314975e-12])
        test_dvdt = acceleration(dummy_r, dummy_m2)
        expected_dvdt = np.around(expected_dvdt, 17)
        test_dvdt = np.around(test_dvdt, 17)
        self.assertEqual(list(expected_dvdt), list(test_dvdt))

    def test_runkut(self):
        ## first timestep
        dummy_rt0 = np.array([1000, 0, 0])
        dummy_vt0 = np.array([0, 100, 0])
        dummy_dt = 0.1
        dummy_m = 100
        # initial accn should be in negative x direction
        vk1 = acceleration(dummy_rt0, dummy_m)
        rk1 = dummy_vt0
        vk2 = acceleration(dummy_rt0 + (rk1 * (dummy_dt/2)), dummy_m)
        rk2 = dummy_vt0 + vk1 * (dummy_dt/2)
        vk3 = acceleration(dummy_rt0 + (rk2 * (dummy_dt/2)), dummy_m)
        rk3 = dummy_vt0 + vk2 * (dummy_dt/2)
        vk4 = acceleration(dummy_rt0 + (rk3 * dummy_dt), dummy_m)
        rk4 = dummy_vt0 + vk3 * dummy_dt

        accn = (1/6) * (vk1 + 2 * vk2 + 2 * vk3 + vk4)
        vel = (1/6) * (rk1 + 2 * rk2 + 2 * rk3 + rk4)

        self.assertLess(accn[0], 0)
        self.assertLess(accn[1], 0)
        self.assertEqual(accn[2], 0)

        self.assertLess(vel[0], 0)
        self.assertGreater(vel[1], 0)
        self.assertEqual(vel[2], 0)


if __name__ == "__main__":
    unittest.main()
