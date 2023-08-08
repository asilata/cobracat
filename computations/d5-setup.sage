R = RootSystem(CartanType(['D',5]))
R.root_lattice().dynkin_diagram()

central_charge = matrix([[-3,1],[-2,1], [-1,1],[7,1],[3,4]])

change_of_basis = matrix([vector(R.ambient_space()(v)) for v in R.root_lattice().simple_roots()])

central_change_euclidean = change_of_basis.inverse() * central_charge
