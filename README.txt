 - �ʱ� �� ����
������ �������� �ʱⰪ�� Multigrid.h ��ܿ��� ���ǵǾ� �ֽ��ϴ�.
1. INITIAL_N : main �Լ��� n�� �����Ǹ�, �ʱ� resolution�� n by n���� �ǹ��մϴ�.
2. ITERATION : �� coarsing ���������� gauss seidel iteration Ƚ���� �ǹ��մϴ�.
3. CYCLE : V-cycle�� �ݺ��� Ƚ���� �ǹ��մϴ�.
4. THRESHOLD : dependency�� ���� �� �� ���� ���Ͽ� threshold ���� ������ �� �ְ� �����Ͽ����ϴ�.
5. USE_THRESHOLD : threshold ���� �̿��Ͽ� ����ϸ� �ð����⵵�� �����մϴ�. �� �κ��� �ذ��ϱ� ���� THRESHOLD ���� ����� ���� ���̶�� �����ϰ� ������ ����ϴ� �˰����� �߰��Ͽ����ϴ�. �� ���� false�� �δ� ���� ��õ�մϴ�.

 - �Լ� ����
1. int coarse(CSR * A, int size, CSR * I, CSR * T)
 matrix A_h�� �̿��Ͽ� dependency �� ����ϰ� �׿� ���� coarsing algorithm�� �����Ͽ�
���������� array d�� c�� ���մϴ�. �̸� �̿��Ͽ� interpolation operator I�� �� transpose form�� T�� ���մϴ�.
������ I_h, T_h�� I^h_2h, I^2h_h �� �ش��մϴ�.
�� �Լ��� coarsing���� ���õ� element���� size�� return�ϰ� �̸� �̿��Ͽ� ���� level�� operator�� matrix�� size�� ���մϴ�.

2. void restrict(CSR * A, double * f, int size, double * u, double * f_2h, CSR * T)
Gauss Seidel�� �Ἥ A*u=f�� Relaxation �� ��, �� residual�� ���ϰ� �̸� restrict �Ͽ� coarse gird�� f_2h�� ���մϴ�.

3. void coarse_A(CSR *T, CSR *A, CSR *I, CSR *A_2h)
Matrix multiplication�� �����մϴ�. ��, T*A*I ���� A_2h�� �Ҵ��մϴ�. 

 - main ����
main�� Multigrid.cpp�� ���� �Ǿ� �ֽ��ϴ�.
1. coarsing�� restrction�� �ݺ��Ͽ� A_16h, f_16h�� ���մϴ�.
2. bottom level���� criterion�� (10^-6���� �����Ͽ����ϴ�.) �̿��� gauss seidel�� ����� �ݺ� Ƚ���� e_16h�� �ٻ簪�� ���մϴ�.
3. Coarse-gird error�� fine grid�� interpolate�ϴ� ������ ���� find-grid������ �ٻ簪�� �����س����ϴ�.
 (���⼭ I_h, ... I_16h �� ���˴ϴ�.)
4. ���������� e_h�� �̿��Ͽ� u_h�� �����ϸ� �ϳ��� cycle�� �����ϴ�.
5. �ռ� ���� operator���� �̿��ؼ� ���� ������ �ٸ��� �ݺ��մϴ�.
6. �־��� �ݺ� Ƚ����ŭ V-cycle�� �����ϸ�, �������� u_h�� �ٻ簪�� ����ϴ�.

 - ��Ÿ
1. Eigen library ������ compile�� ���� �ʴ� ���, ��θ� �����Ͽ� �־�� �մϴ�.
visual C++ 2010 express ��� ������ ���� ������� �ذ��Ͽ����ϴ�.
���� project(Multigrid) ��Ŭ�� - properties - VC++ Directories - Include Directories -> eigen ���� ��ο� �߰�.

2. Coarsing�� detail�� ������ �������� �ʾ�����, Multigrid.h�� ���ǵǾ� �ִ� build_d, build_c ���� �Լ��� �̿��մϴ�.
Coarsing �������� �˰����� �ð����⵵�� ���̱� ���ؼ� index��� �迭�� ����Ͽ� �� element�鿡 dependent�� element���� index�� �����Ͽ����ϴ�.