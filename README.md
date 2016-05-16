##Algebraic multigrid method

�־��� f(x,y) = 2[(1-6x<sup>2</sup>)y<sup>2</sup>(1-y<sup>2</sup>) + (1-6y<sup>2</sup>)x<sup>2</sup>(1-x<sup>2</sup>)]�� ����,  
�̺� ������ -u<sub>xx</sub>-u<sub>yy</sub> = f(x,y)�� �ַ�� u�� ���ϴ� ���� ������Ʈ�� ��ǥ�̴�.  
(domain�� [0,1]x[0,1]�̸�, boundary condition�� u=0�̴�.)  

�̺� �������� Ǯ�� Poisson's equation���κ��� Au=h�� ��´�.  
���⼭ u�� n<sup>2</sup> by 1 vector�̰�, A��  n<sup>2</sup> by  n<sup>2</sup> matrix�� �ȴ�.  
�� ������Ʈ������ �־��� Possion's equation�� Ǯ�� ���� Algebraic MultiGrid(AMG) method�� �����Ͽ� A�� inverse�� ���Ѵ�.  

�� ������Ʈ�� 'Microsoft Visual C++ 2010 Express'�� ����Ͽ�����, compile�� execution file�� 'results' folder���� Ȯ���� �� �ִ�.  
execution file�� �����ϸ�, solution vector u�� 'output.txt'�� �����Ѵ�.  
matlab function *'drawResult'*�� ����Ͽ� ����� Ȯ���� �� �ִ�.  


#### �ʱ� ��
������ �������� �ʱⰪ�� Multigrid.h ��ܿ��� ���ǵǾ� �ֽ��ϴ�.
<ol>
<li> INITIAL_N : main �Լ��� n�� �����Ǹ�, �ʱ� resolution�� n by n���� �ǹ��մϴ�.</li>
<li> ITERATION : �� coarsing ���������� gauss seidel iteration Ƚ���� �ǹ��մϴ�.</li>
<li> CYCLE : V-cycle�� �ݺ��� Ƚ���� �ǹ��մϴ�.</li>
<li> THRESHOLD : dependency�� ���� �� �� ���� ���Ͽ� threshold ���� ������ �� �ְ� �����Ͽ����ϴ�.</li>
<li> USE_THRESHOLD : threshold ���� �̿��Ͽ� ����ϸ� �ð����⵵�� �����մϴ�. �� �κ��� �ذ��ϱ� ���� THRESHOLD ���� ����� ���� ���̶�� �����ϰ� ������ ����ϴ� �˰����� �߰��Ͽ����ϴ�. �� ���� false�� �δ� ���� ��õ�մϴ�.</li>
</ol>

#### �Լ�
<ol>
<li> int coarse(CSR * A, int size, CSR * I, CSR * T)</li>
<p>matrix A_h�� �̿��Ͽ� dependency �� ����ϰ� �׿� ���� coarsing algorithm�� �����Ͽ� ���������� array d�� c�� ���մϴ�. �̸� �̿��Ͽ� interpolation operator I�� �� transpose form�� T�� ���մϴ�.
������ I_h, T_h�� I^h_2h, I^2h_h �� �ش��մϴ�. �� �Լ��� coarsing���� ���õ� element���� size�� return�ϰ� �̸� �̿��Ͽ� ���� level�� operator�� matrix�� size�� ���մϴ�.</p>

<li> void restrict(CSR * A, double * f, int size, double * u, double * f_2h, CSR * T)</li>
<p>Gauss Seidel�� �Ἥ A*u=f�� Relaxation �� ��, �� residual�� ���ϰ� �̸� restrict �Ͽ� coarse gird�� f_2h�� ���մϴ�.</p>

<li> void coarse_A(CSR *T, CSR *A, CSR *I, CSR *A_2h)</li>
<p>Matrix multiplication�� �����մϴ�. ��, T*A*I ���� A_2h�� �Ҵ��մϴ�. </p>
</ol>

#### main
main�� Multigrid.cpp�� ���� �Ǿ� �ֽ��ϴ�.
<ol>
<li> coarsing�� restrction�� �ݺ��Ͽ� A_16h, f_16h�� ���մϴ�. </li>
<li> bottom level���� criterion�� (10^-6���� �����Ͽ����ϴ�.) �̿��� gauss seidel�� ����� �ݺ� Ƚ���� e_16h�� �ٻ簪�� ���մϴ�.</li>
<li> Coarse-gird error�� fine grid�� interpolate�ϴ� ������ ���� find-grid������ �ٻ簪�� �����س����ϴ�.
 (���⼭ I_h, ... I_16h �� ���˴ϴ�.)</li>
<li> ���������� e_h�� �̿��Ͽ� u_h�� �����ϸ� �ϳ��� cycle�� �����ϴ�.</li>
<li> �ռ� ���� operator���� �̿��ؼ� ���� ������ �ٸ��� �ݺ��մϴ�.</li>
<li> �־��� �ݺ� Ƚ����ŭ V-cycle�� �����ϸ�, �������� u_h�� �ٻ簪�� ����ϴ�.</li>
</ol>

#### ��Ÿ
<ul>
<li>Eigen library ������ compile�� ���� �ʴ� ���, ��θ� �����Ͽ� �־�� �մϴ�.</li>
<p>visual C++ 2010 express ��� ������ ���� ������� �ذ��Ͽ����ϴ�.</p>
<p>���� project(Multigrid) ��Ŭ�� - properties - VC++ Directories - Include Directories -> eigen ���� ��ο� �߰�.</p>

<li>Coarsing�� detail�� ������ �������� �ʾ�����, Multigrid.h�� ���ǵǾ� �ִ� build_d, build_c ���� �Լ��� �̿��մϴ�.</li>
<p>Coarsing �������� �˰����� �ð����⵵�� ���̱� ���ؼ� index��� �迭�� ����Ͽ� �� element�鿡 dependent�� element���� index�� �����Ͽ����ϴ�.</p>
</ul>
