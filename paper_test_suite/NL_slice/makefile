# NL_slice tests run for the paper

convergence_order_1:
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=1 --ncells_1d=50
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=1 --ncells_1d=60
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=1 --ncells_1d=70
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=1 --ncells_1d=80
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=1 --ncells_1d=90
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=1 --ncells_1d=100
    python NL_slice.py --configuration=convergence --order=1 --ncells_1d=50
    python NL_slice.py --configuration=convergence --order=1 --ncells_1d=60
    python NL_slice.py --configuration=convergence --order=1 --ncells_1d=70
    python NL_slice.py --configuration=convergence --order=1 --ncells_1d=80
    python NL_slice.py --configuration=convergence --order=1 --ncells_1d=90
    python NL_slice.py --configuration=convergence --order=1 --ncells_1d=100

convergence_order_0:
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=0 --ncells_1d=100
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=0 --ncells_1d=120
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=0 --ncells_1d=140
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=0 --ncells_1d=160
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=0 --ncells_1d=180
    python NL_slice.py --conservative_transport=True --configuration=convergence --order=0 --ncells_1d=200
    python NL_slice.py --configuration=convergence --order=0 --ncells_1d=100
    python NL_slice.py --configuration=convergence --order=0 --ncells_1d=120
    python NL_slice.py --configuration=convergence --order=0 --ncells_1d=140
    python NL_slice.py --configuration=convergence --order=0 --ncells_1d=160
    python NL_slice.py --configuration=convergence --order=0 --ncells_1d=180
    python NL_slice.py --configuration=convergence --order=0 --ncells_1d=200

convergence_plotting:
    python plot_NL_slice_error_and_Td.py

consistency:
    python NL_slice.py --configuration=consistency --order=1 --ncells_1d=100
    python NL_slice.py --conservative_transport=True --configuration=consistency --order=1 --ncells_1d=100
    python NL_slice.py --configuration=consistency --order=0 --ncells_1d=200
    python NL_slice.py --conservative_transport=True --configuration=consistency --order=0 --ncells_1d=200

consistency_plotting:
    python plot_NL_slice_consistency_mx_sse.py