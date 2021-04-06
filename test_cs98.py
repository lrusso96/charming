from charm.toolbox.eccurve import prime192v1
from cs98 import CS98, InputError

import unittest


class TestCS98(unittest.TestCase):

    def test_basic(self):
        curve = prime192v1
        secparam = 128

        cs = CS98(curve)
        pk, sk = cs.keygen(secparam)
        msg = b"Message of length 20"
        c = cs.encrypt(pk, msg)
        m = cs.decrypt(pk, sk, c)
        self.assertEqual(m, msg, "Should be equal")

    def test_short_msg(self):
        curve = prime192v1
        secparam = 128

        cs = CS98(curve)
        pk, _ = cs.keygen(secparam)

        msg = b"Short msg"
        self.assertRaises(InputError, cs.encrypt, pk, msg)

    def test_long_msg(self):
        curve = prime192v1
        secparam = 128

        cs = CS98(curve)
        pk, _ = cs.keygen(secparam)

        msg = b"Too long message should not work!"
        self.assertRaises(InputError, cs.encrypt, pk, msg)


if __name__ == '__main__':
    unittest.main()
