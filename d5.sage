load("stability-conditions.sage")

cc = CentralCharge("D5", (-9 + 4*I, -8 + 5*I, 9 + 5*I, -3 + 4*I, 9 + 2*I))
sc = StandardStabilityCondition(cc)
alpha = cc.simple_roots

# 9 element cell
