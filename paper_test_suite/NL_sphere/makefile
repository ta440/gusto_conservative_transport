# Documenting the simulations run for the NL_sphere tests for the paper.

convergence_order_1:
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=8
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=16
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=24
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=1 --ncells_per_edge=32

	python NL_sphere_divergent.py --configuration=convergence --order=1 --ncells_per_edge=8
	python NL_sphere_divergent.py --configuration=convergence --order=1 --ncells_per_edge=16
	python NL_sphere_divergent.py --configuration=convergence --order=1 --ncells_per_edge=24
	python NL_sphere_divergent.py --configuration=convergence --order=1 --ncells_per_edge=32

convergence_order_0:
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=16
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=24
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=36
	python NL_sphere_divergent.py --conservative_transport=True --configuration=convergence --order=0 --ncells_per_edge=48

	python NL_sphere_divergent.py --configuration=convergence --order=0 --ncells_per_edge=16
	python NL_sphere_divergent.py --configuration=convergence --order=0 --ncells_per_edge=24
	python NL_sphere_divergent.py --configuration=convergence --order=0 --ncells_per_edge=36
	python NL_sphere_divergent.py --configuration=convergence --order=0 --ncells_per_edge=48

convergence plotting:
	python plot_NL_sphere_error_and_Td.py

consistency:
	python NL_sphere_divergent.py --configuration=consistency --order=1 --ncells_per_edge=32
	python NL_sphere_divergent.py --conservative_transport=True --configuration=consistency --order=1 --ncells_per_edge=32
	python NL_sphere_divergent.py --configuration=consistency --order=0 --ncells_per_edge=48
	python NL_sphere_divergent.py --conservative_transport=True --configuration=consistency --order=0 --ncells_per_edge=48

consistency plotting:
	python plot_NL_sphere_consistency_mx_sse.py
