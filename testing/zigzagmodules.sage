load("../zigzagalgebra.sage")
load("../zigzagmodules.sage")
R = make_test(a3graph)['Z']
R.inject_variables()

P = ZigZagModule(R, R.idempotents()[0])

