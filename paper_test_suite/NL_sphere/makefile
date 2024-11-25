# Documenting the simulations run for the paper.
# Compare spatial resolutions of ncells per edge of 8,16,32,
# which correspond to reference levels of 3,4,5. 

convergence_order_1:
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=8
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=16
	python NL_sphere.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=32

	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=8
	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=16
	python NL_sphere.py --configuration=convergence --order=1 --ncells_per_edge=32

convergence_order_0:
