 - 초기 값 설정
다음의 변수들의 초기값은 Multigrid.h 상단에서 정의되어 있습니다.
1. INITIAL_N : main 함수의 n에 대응되며, 초기 resolution이 n by n임을 의미합니다.
2. ITERATION : 각 coarsing 과정에서의 gauss seidel iteration 횟수를 의미합니다.
3. CYCLE : V-cycle을 반복할 횟수를 의미합니다.
4. THRESHOLD : dependency를 구할 때 이 값을 통하여 threshold 값을 지정할 수 있게 설정하였습니다.
5. USE_THRESHOLD : threshold 값을 이용하여 계산하면 시간복잡도가 증가합니다. 이 부분을 해결하기 위해 THRESHOLD 값이 충분히 작은 값이라고 가정하고 빠르게 계산하는 알고리즘을 추가하였습니다. 이 값을 false로 두는 것을 추천합니다.

 - 함수 설명
1. int coarse(CSR * A, int size, CSR * I, CSR * T)
 matrix A_h를 이용하여 dependency 를 계산하고 그에 따라 coarsing algorithm을 적용하여
내부적으로 array d와 c를 구합니다. 이를 이용하여 interpolation operator I와 그 transpose form인 T를 구합니다.
각각의 I_h, T_h는 I^h_2h, I^2h_h 에 해당합니다.
이 함수는 coarsing으로 선택된 element들의 size를 return하고 이를 이용하여 하위 level의 operator와 matrix의 size를 정합니다.

2. void restrict(CSR * A, double * f, int size, double * u, double * f_2h, CSR * T)
Gauss Seidel을 써서 A*u=f를 Relaxation 한 뒤, 그 residual을 구하고 이를 restrict 하여 coarse gird의 f_2h를 구합니다.

3. void coarse_A(CSR *T, CSR *A, CSR *I, CSR *A_2h)
Matrix multiplication을 수행합니다. 즉, T*A*I 값을 A_2h에 할당합니다. 

 - main 설명
main은 Multigrid.cpp에 정의 되어 있습니다.
1. coarsing과 restrction을 반복하여 A_16h, f_16h를 구합니다.
2. bottom level에서 criterion을 (10^-6으로 설정하였습니다.) 이용한 gauss seidel로 충분한 반복 횟수로 e_16h의 근사값을 구합니다.
3. Coarse-gird error를 fine grid로 interpolate하는 과정을 통해 find-grid에서의 근사값을 수정해나갑니다.
 (여기서 I_h, ... I_16h 가 사용됩니다.)
4. 최종적으로 e_h를 이용하여 u_h를 보완하면 하나의 cycle이 끝납니다.
5. 앞서 구한 operator들을 이용해서 위의 과정을 바르게 반복합니다.
6. 주어진 반복 횟수만큼 V-cycle을 시행하면, 최종적인 u_h의 근사값을 얻습니다.

 - 기타
1. Eigen library 문제로 compile이 되지 않는 경우, 경로를 설정하여 주어야 합니다.
visual C++ 2010 express 경우 다음과 같은 방법으로 해결하였습니다.
좌측 project(Multigrid) 우클릭 - properties - VC++ Directories - Include Directories -> eigen 폴더 경로에 추가.

2. Coarsing의 detail한 과정은 설명하지 않았으나, Multigrid.h에 정의되어 있는 build_d, build_c 등의 함수를 이용합니다.
Coarsing 과정에서 알고리즘의 시간복잡도를 줄이기 위해서 index라는 배열을 사용하여 각 element들에 dependent한 element들의 index를 저장하였습니다.