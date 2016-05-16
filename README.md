##Algebraic multigrid method

주어진 f(x,y) = 2[(1-6x<sup>2</sup>)y<sup>2</sup>(1-y<sup>2</sup>) + (1-6y<sup>2</sup>)x<sup>2</sup>(1-x<sup>2</sup>)]에 대해,  
미분 방정식 -u<sub>xx</sub>-u<sub>yy</sub> = f(x,y)의 솔루션 u를 구하는 것이 프로젝트의 목표이다.  
(domain은 [0,1]x[0,1]이며, boundary condition은 u=0이다.)  

미분 방정식을 풀면 Poisson's equation으로부터 Au=h를 얻는다.  
여기서 u는 n<sup>2</sup> by 1 vector이고, A는  n<sup>2</sup> by  n<sup>2</sup> matrix가 된다.  
본 프로젝트에서는 주어진 Possion's equation을 풀기 위해 Algebraic MultiGrid(AMG) method를 적용하여 A의 inverse를 구한다.  

본 프로젝트는 'Microsoft Visual C++ 2010 Express'을 사용하였으며, compile된 execution file은 'results' folder에서 확인할 수 있다.  
execution file을 실행하면, solution vector u를 'output.txt'로 저장한다.  
matlab function *'drawResult'*를 사용하여 결과를 확인할 수 있다.  


#### 초기 값
다음의 변수들의 초기값은 Multigrid.h 상단에서 정의되어 있습니다.
<ol>
<li> INITIAL_N : main 함수의 n에 대응되며, 초기 resolution이 n by n임을 의미합니다.</li>
<li> ITERATION : 각 coarsing 과정에서의 gauss seidel iteration 횟수를 의미합니다.</li>
<li> CYCLE : V-cycle을 반복할 횟수를 의미합니다.</li>
<li> THRESHOLD : dependency를 구할 때 이 값을 통하여 threshold 값을 지정할 수 있게 설정하였습니다.</li>
<li> USE_THRESHOLD : threshold 값을 이용하여 계산하면 시간복잡도가 증가합니다. 이 부분을 해결하기 위해 THRESHOLD 값이 충분히 작은 값이라고 가정하고 빠르게 계산하는 알고리즘을 추가하였습니다. 이 값을 false로 두는 것을 추천합니다.</li>
</ol>

#### 함수
<ol>
<li> int coarse(CSR * A, int size, CSR * I, CSR * T)</li>
<p>matrix A_h를 이용하여 dependency 를 계산하고 그에 따라 coarsing algorithm을 적용하여 내부적으로 array d와 c를 구합니다. 이를 이용하여 interpolation operator I와 그 transpose form인 T를 구합니다.
각각의 I_h, T_h는 I^h_2h, I^2h_h 에 해당합니다. 이 함수는 coarsing으로 선택된 element들의 size를 return하고 이를 이용하여 하위 level의 operator와 matrix의 size를 정합니다.</p>

<li> void restrict(CSR * A, double * f, int size, double * u, double * f_2h, CSR * T)</li>
<p>Gauss Seidel을 써서 A*u=f를 Relaxation 한 뒤, 그 residual을 구하고 이를 restrict 하여 coarse gird의 f_2h를 구합니다.</p>

<li> void coarse_A(CSR *T, CSR *A, CSR *I, CSR *A_2h)</li>
<p>Matrix multiplication을 수행합니다. 즉, T*A*I 값을 A_2h에 할당합니다. </p>
</ol>

#### main
main은 Multigrid.cpp에 정의 되어 있습니다.
<ol>
<li> coarsing과 restrction을 반복하여 A_16h, f_16h를 구합니다. </li>
<li> bottom level에서 criterion을 (10^-6으로 설정하였습니다.) 이용한 gauss seidel로 충분한 반복 횟수로 e_16h의 근사값을 구합니다.</li>
<li> Coarse-gird error를 fine grid로 interpolate하는 과정을 통해 find-grid에서의 근사값을 수정해나갑니다.
 (여기서 I_h, ... I_16h 가 사용됩니다.)</li>
<li> 최종적으로 e_h를 이용하여 u_h를 보완하면 하나의 cycle이 끝납니다.</li>
<li> 앞서 구한 operator들을 이용해서 위의 과정을 바르게 반복합니다.</li>
<li> 주어진 반복 횟수만큼 V-cycle을 시행하면, 최종적인 u_h의 근사값을 얻습니다.</li>
</ol>

#### 기타
<ul>
<li>Eigen library 문제로 compile이 되지 않는 경우, 경로를 설정하여 주어야 합니다.</li>
<p>visual C++ 2010 express 경우 다음과 같은 방법으로 해결하였습니다.</p>
<p>좌측 project(Multigrid) 우클릭 - properties - VC++ Directories - Include Directories -> eigen 폴더 경로에 추가.</p>

<li>Coarsing의 detail한 과정은 설명하지 않았으나, Multigrid.h에 정의되어 있는 build_d, build_c 등의 함수를 이용합니다.</li>
<p>Coarsing 과정에서 알고리즘의 시간복잡도를 줄이기 위해서 index라는 배열을 사용하여 각 element들에 dependent한 element들의 index를 저장하였습니다.</p>
</ul>
