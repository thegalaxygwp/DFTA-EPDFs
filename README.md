# DFTA
This fortran project proposes a method that constructs linear equations system through discrete Fourier series and solves the elliptic PDEs by inverse matrix of its coefficient. 
In the algorithm design, we make full use of boundary conditions and research mode numbers to simplify the solution process and reduce computational workload. Such as sampling less grid points under a linear sigle mode. 
It is particularly applicable when solving repeatedly while the differential coefficient is specific and invariant. Application shows that this method can provide extremly high accuracy under sufficient series decomposition. 
It use function collection like intel-MKL, LAPACK library and OpenMP.
