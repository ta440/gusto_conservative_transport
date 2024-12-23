# Documenting the simulations run for the paper.
# Compare spatial resolutions of ncells per edge of 6,12,24,48 
# this is the number of cells on edge of a cubed sphere face

convergence_order_1:
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=6
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=12
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=24
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=48

	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=6
	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=12
	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=24
	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=48

convergence_order_0:
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=6
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=12
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=24
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=48

	python NL_sphere.py --configuration=convergence --order=0 --ncells_per_edge=6
	python NL_sphere.py --configuration=convergence --order=0 --ncells_per_edge=12
	python NL_sphere.py --configuration=convergence --order=0 --ncells_per_edge=24
	python NL_sphere.py --configuration=convergence --order=0 --ncells_per_edge=48
