#include<iostream>
#include<vector>
#include<string>

using namespace std;

//const double current = 10.0;
int main()
{
	cout << "���ڿ�ʼ����ά����Ч��·" << endl;

	//�����ж��ٸ��ڵ�
	cout << "������ڵ����" << endl;
	int nodeNums = 0;
	cin >> nodeNums;

	cout << "������Ҫ���Ľڵ�" << endl;
	int ansNode1, ansNode2;
	cin >> ansNode1 >> ansNode2;

	//����˵�
	cout << "This is the meun" << endl;
	cout << "e--->�������             u--->����ڵ��ĵ�ѹ" << endl;
	cout << "i--->����ڵ��ĵ���      r--->����ڵ�֮��ĵ�ѹ" << endl;
	cout << "o--->�����ܿ�Դģʽ" << endl;

	string c = "";
	cin >> c;

	//��������A�����絼ϵ����������ȫ����ɺ�ѡȡ���˵�һ�к͵�һ�е�����Ԫ������µľ���
	vector<vector<double>> A(nodeNums);
	for (int i = 0;i < nodeNums;i++)
		for (int j = 0;j < nodeNums;j++)
			(A[i]).push_back(0.0);

	//������������B
	vector<double> B(nodeNums, 0.0);

	while (c != "e")
	{
		if (c == "r")
		{
			cout << "���������ڵ��ź͵���ֵ" << endl;
			string flag = "y";
			while (flag == "y")
			{
				int node1, node2;
				double g;
				cin >> node1 >> node2 >> g;

				g = 1.0 / g;
				A[node1][node1] += g;
				A[node1][node2] -= g;
				A[node2][node1] -= g;
				A[node2][node2] += g;

				cout << "�Ƿ�Ҫ������裿��y/n)" << endl;
				cin >> flag;
			}
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
		if (c == "u")
		{
			//��ȡ��ѹԴ��ֵ
			cout << "�������ѹԴ�ڵ��ź͵�ѹֵ(��һ��Ϊ�ߵ�ѹ�ڵ㣬�ڶ���Ϊ�͵�ѹ)" << endl;
			string flag = "y";
			while (flag == "y")
			{
				int p, n;
				double u;
				cin >> p >> n >> u;

				int m = A.size();

				A.push_back(vector<double>(A.size(), 0.0));
				for (int i = 0;i < A.size();i++)
					A[i].push_back(0.0);
				B.push_back(0.0);

				A[p][m] += 1;
				A[n][m] -= 1;
				A[m][p] += 1;
				A[m][n] -= 1;
				B[m] += u;

				cout << "�Ƿ�Ҫ�����ѹ����y/n)" << endl;
				cin >> flag;
			}
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
		if (c == "i")
		{
			cout << "����������ڵ��ź͵���ֵ" << endl;

			string flag = "y";
			while (flag == "y")
			{
				int node1, node2;
				double i;
				cin >> node1 >> node2 >> i;

				B[node1] -= i;
				B[node2] += i;

				cout << "�Ƿ�Ҫ�����������y/n)" << endl;
				cin >> flag;
			}
			//B[ansNode1] -= current;
			//B[ansNode2] += current;
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
		if (c == "o")
		{
			cout << "Ŀǰ�������ܿ�Դ�ɹ�ѡ��(�����ǲ˵�)" << endl;
			cout << "1.��ѹ���Ƶĵ���Դ()" << endl
				<< "2.��ѹ���Ƶĵ�ѹԴ()" << endl
				<< "3.�������Ƶĵ���Դ()" << endl;

			cout << "���������ѡ��" << endl;

			string choice = "";
			cin >> choice;
			//ϸ��ѡ��
			if (choice == "1")
			{
				cout << "�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������" << endl;
				string flag = "y";
				while (flag == "y")
				{
					int p, n, pc, nc;
					double g;
					cin >> p >> n >> pc >> nc >> g;

					A[p][pc] += g;
					A[p][nc] -= g;
					A[n][pc] -= g;
					A[n][nc] += g;

					cout << "�Ƿ�Ҫ���루y/n��" << endl;
					cin >> flag;
				}
			}
			if (choice == "2")
			{
				cout << "�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������" << endl;
				string flag = "y";
				while (flag == "y")
				{
					int p, n, pc, nc;
					double a;
					cin >> p >> n >> pc >> nc >> a;

					int m = A.size();
					A.push_back(vector<double>(A.size(), 0.0));
					for (int i = 0;i < A.size();i++)
						A[i].push_back(0.0);
					B.push_back(0.0);
					A[p][m] += 1;
					A[n][m] -= 1;
					A[m][p] += 1;
					A[m][n] -= 1;
					A[m][pc] -= a;
					A[m][nc] += a;

					cout << "�Ƿ�Ҫ���루y/n��" << endl;
					cin >> flag;
				}
			}
			if (choice == "3")
			{
				cout << "�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������" << endl;
				string flag = "y";
				while (flag == "y")
				{
					int p, n, pc, nc;
					double b;
					cin >> p >> n >> pc >> nc >> b;

					int m = A.size();
					A.push_back(vector<double>(A.size(), 0.0));
					for (int i = 0;i < A.size();i++)
						A[i].push_back(0.0);

					B.push_back(0.0);
					A[p][m] += b;
					A[n][m] -= b;

					cout << "�Ƿ�Ҫ���루y/n��" << endl;
					cin >> flag;
				}
			}
			/*
			if(choice=="4")
			{
				cout<<"�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������"<<endl;
				string flag="y";
				while(flag=="y")
				{
					int p, n, pc, nc;
					double r;
					cin >> p >> n >> pc >> nc >> r;

					A.push_back(vector<double>(A.size(), 0.0));
					for (int i = 0;i < A.size();i++)
						A[i].push_back(0.0);

					int m = A.size();
					int k = 0;
					A[p][m] += 1;
					A[n][m] -= 1;
					A[m][p] += 1;
					A[m][n] -= 1;
					A[m][k] -= r;
				}
			}
			*/
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
	}

	//��������������
	nodeNums = A.size();
	vector<vector<double>> AN(nodeNums - 1);
	for (int i = 0;i < nodeNums - 1;i++)
		for (int j = 0;j < nodeNums - 1;j++)
			AN[i].push_back(A[i + 1][j + 1]);
	vector<double> BN;
	for (int i = 0;i < nodeNums - 1;i++)
		BN.push_back(B[i + 1]);

	//�������Ľ��
	//LU����ֽ�
	vector<vector<double>> l(nodeNums - 1);
	for (int i = 0;i < nodeNums - 1;i++)
		for (int j = 0;j < nodeNums - 1;j++)
			(l[i]).push_back(0.0);
	vector<vector<double>> u(nodeNums - 1);
	for (int i = 0;i < nodeNums - 1;i++)
		for (int j = 0;j < nodeNums - 1;j++)
			(u[i]).push_back(0.0);

	int i, r, k;

	const int n = nodeNums - 1;

	//����U�ĵ�һ�еĸ�ֵ
	for (i = 0; i < n; i++)
		u[0][i] = AN[0][i];

	//����L�ĵ�һ�еĸ�ֵ
	for (i = 1; i < n; i++)
		l[i][0] = AN[i][0] / u[0][0];

	//����U��ʣ�µ�������L��ʣ�µ�����
	for (r = 1; r < n; r++)
	{
		for (i = r; i < n; i++)
		{
			double sum1 = 0;
			for (k = 0; k < r; k++)
				sum1 += l[r][k] * u[k][i];
			u[r][i] = AN[r][i] - sum1;
		}

		if (r != n)
			for (i = r + 1;i < n;i++)
			{
				double sum2 = 0;
				for (k = 0; k < r; k++)
					sum2 += l[i][k] * u[k][r];
				l[i][r] = (AN[i][r] - sum2) / u[r][r];
			}
	}

	vector<double> y(n, 0.0);
	y[0] = BN[0];
	for (i = 1; i < n; i++)
	{
		double sum3 = 0;
		for (k = 0; k < i; k++)
			sum3 += l[i][k] * y[k];
		y[i] = BN[i] - sum3;
	}
	vector<double> x(n, 0.0);
	x[n - 1] = y[n - 1] / u[n - 1][n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		double sum4 = 0;
		for (k = i + 1; k < n; k++)
			sum4 += u[i][k] * x[k];
		x[i] = (y[i] - sum4) / u[i][i];
	}
	vector<double> u_v(n+1,0.0);
	u_v[0]=0.0;
	for(int i=1;i<n+1;i++)
		u_v[i]=x[i-1];
	double ansu = u_v[ansNode1] - u_v[ansNode2];


	//��������Ч����

	//�����ж��ٸ��ڵ�
	cout << "������ڵ����" << endl;
	nodeNums = 0;
	cin >> nodeNums;
	double input_i = 1.0;
	/*
	cout << "������Ҫ���Ľڵ�" << endl;
	int ansNode1, ansNode2;
	cin >> ansNode1 >> ansNode2;
	*/
	//����˵�
	cout << "This is the meun" << endl;
	cout << "e--->�������             u--->����ڵ��ĵ�ѹ" << endl;
	cout << "i--->����ڵ��ĵ���      r--->����ڵ�֮��ĵ�ѹ" << endl;
	cout << "o--->�����ܿ�Դģʽ" << endl;

	c = "";
	cin >> c;

	//��������A�����絼ϵ����������ȫ����ɺ�ѡȡ���˵�һ�к͵�һ�е�����Ԫ������µľ���
	vector<vector<double>> AA(nodeNums);
	for (int i = 0;i < nodeNums;i++)
		for (int j = 0;j < nodeNums;j++)
			(AA[i]).push_back(0.0);

	//������������B
	vector<double> BB(nodeNums, 0.0);

	while (c != "e")
	{
		if (c == "r")
		{
			cout << "���������ڵ��ź͵���ֵ" << endl;
			string flag = "y";
			while (flag == "y")
			{
				int node1, node2;
				double g;
				cin >> node1 >> node2 >> g;

				g = 1.0 / g;
				AA[node1][node1] += g;
				AA[node1][node2] -= g;
				AA[node2][node1] -= g;
				AA[node2][node2] += g;

				cout << "�Ƿ�Ҫ������裿��y/n)" << endl;
				cin >> flag;
			}
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
		if (c == "u")
		{
			//��ȡ��ѹԴ��ֵ
			cout << "�������ѹԴ�ڵ��ź͵�ѹֵ(��һ��Ϊ�ߵ�ѹ�ڵ㣬�ڶ���Ϊ�͵�ѹ)" << endl;
			string flag = "y";
			while (flag == "y")
			{
				int p, n;
				double u;
				cin >> p >> n >> u;

				int m = AA.size();

				AA.push_back(vector<double>(AA.size(), 0.0));
				for (int i = 0;i < AA.size();i++)
					AA[i].push_back(0.0);
				BB.push_back(0.0);

				AA[p][m] += 1;
				AA[n][m] -= 1;
				AA[m][p] += 1;
				AA[m][n] -= 1;
				BB[m] += u;

				cout << "�Ƿ�Ҫ�����ѹ����y/n)" << endl;
				cin >> flag;
			}
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
		if (c == "i")
		{
			cout << "����������ڵ��ź͵���ֵ" << endl;

			string flag = "y";
			while (flag == "y")
			{
				int node1, node2;
				double i;
				cin >> node1 >> node2 >> i;
				input_i = i;
				BB[node1] -= i;
				BB[node2] += i;

				cout << "�Ƿ�Ҫ�����������y/n)" << endl;
				cin >> flag;
			}
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
		if (c == "o")
		{
			cout << "Ŀǰ�������ܿ�Դ�ɹ�ѡ��(�����ǲ˵�)" << endl;
			cout << "1.��ѹ���Ƶĵ���Դ()" << endl
				<< "2.��ѹ���Ƶĵ�ѹԴ()" << endl
				<< "3.�������Ƶĵ���Դ()" << endl;

			cout << "���������ѡ��" << endl;

			string choice = "";
			cin >> choice;
			//ϸ��ѡ��
			if (choice == "1")
			{
				cout << "�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������" << endl;
				string flag = "y";
				while (flag == "y")
				{
					int p, n, pc, nc;
					double g;
					cin >> p >> n >> pc >> nc >> g;

					AA[p][pc] += g;
					AA[p][nc] -= g;
					AA[n][pc] -= g;
					AA[n][nc] += g;

					cout << "�Ƿ�Ҫ���루y/n��" << endl;
					cin >> flag;
				}
			}
			if (choice == "2")
			{
				cout << "�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������" << endl;
				string flag = "y";
				while (flag == "y")
				{
					int p, n, pc, nc;
					double a;
					cin >> p >> n >> pc >> nc >> a;

					int m = AA.size();
					AA.push_back(vector<double>(AA.size(), 0.0));
					for (int i = 0;i < AA.size();i++)
						AA[i].push_back(0.0);
					BB.push_back(0.0);
					AA[p][m] += 1;
					AA[n][m] -= 1;
					AA[m][p] += 1;
					AA[m][n] -= 1;
					AA[m][pc] -= a;
					AA[m][nc] += a;

					cout << "�Ƿ�Ҫ���루y/n��" << endl;
					cin >> flag;
				}
			}
			if (choice == "3")
			{
				cout << "�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������" << endl;
				string flag = "y";
				while (flag == "y")
				{
					int p, n, pc, nc;
					double b;
					cin >> p >> n >> pc >> nc >> b;

					int m = AA.size();
					AA.push_back(vector<double>(AA.size(), 0.0));
					for (int i = 0;i < AA.size();i++)
						AA[i].push_back(0.0);

					BB.push_back(0.0);
					AA[p][m] += b;
					AA[n][m] -= b;

					cout << "�Ƿ�Ҫ���루y/n��" << endl;
					cin >> flag;
				}
			}
			/*
			if(choice=="4")
			{
				cout<<"�밴�տ��ƽڵ㣬�ܿؽڵ㣬��������˳���������"<<endl;
				string flag="y";
				while(flag=="y")
				{
					int p, n, pc, nc;
					double r;
					cin >> p >> n >> pc >> nc >> r;

					A.push_back(vector<double>(A.size(), 0.0));
					for (int i = 0;i < A.size();i++)
						A[i].push_back(0.0);

					int m = A.size();
					int k = 0;
					A[p][m] += 1;
					A[n][m] -= 1;
					A[m][p] += 1;
					A[m][n] -= 1;
					A[m][k] -= r;
				}
			}
			*/
			cout << "������ѡ�����e����" << endl;
			cin >> c;
		}
	}

	//��������������
	nodeNums = AA.size();
	vector<vector<double>> AAN(nodeNums - 1);
	for (int i = 0;i < nodeNums - 1;i++)
		for (int j = 0;j < nodeNums - 1;j++)
			AAN[i].push_back(AA[i + 1][j + 1]);
	vector<double> BBN;
	for (int i = 0;i < nodeNums - 1;i++)
		BBN.push_back(BB[i + 1]);

	//�������Ľ��
	//LU����ֽ�
	vector<vector<double>> ll(nodeNums - 1);
	for (int i = 0;i < nodeNums - 1;i++)
		for (int j = 0;j < nodeNums - 1;j++)
			(ll[i]).push_back(0.0);
	vector<vector<double>> uu(nodeNums - 1);
	for (int i = 0;i < nodeNums - 1;i++)
		for (int j = 0;j < nodeNums - 1;j++)
			(uu[i]).push_back(0.0);

	int ii, rr, kk;

	const int nn = nodeNums - 1;

	//����U�ĵ�һ�еĸ�ֵ
	for (ii = 0; ii < nn; ii++)
		uu[0][ii] = AAN[0][ii];

	//����L�ĵ�һ�еĸ�ֵ
	for (ii = 1; ii < nn; ii++)
		ll[ii][0] = AAN[ii][0] / uu[0][0];

	//����U��ʣ�µ�������L��ʣ�µ�����
	for (rr = 1; rr < nn; rr++)
	{
		for (ii = rr; ii < nn; ii++)
		{
			double sum1 = 0;
			for (kk = 0; kk < rr; kk++)
				sum1 += ll[rr][kk] * uu[kk][ii];
			uu[rr][ii] = AAN[rr][ii] - sum1;
		}

		if (rr != nn)
			for (ii = rr + 1;ii < nn;ii++)
			{
				double sum2 = 0;
				for (kk = 0; kk < rr; kk++)
					sum2 += ll[ii][kk] * uu[kk][rr];
				ll[ii][rr] = (AAN[ii][rr] - sum2) / uu[rr][rr];
			}
	}

	vector<double> yy(n, 0.0);
	yy[0] = BBN[0];
	for (ii = 1; ii < nn; ii++)
	{
		double sum3 = 0;
		for (kk = 0; kk < ii; kk++)
			sum3 += ll[ii][kk] * yy[kk];
		yy[ii] = BBN[ii] - sum3;
	}
	vector<double> xx(nn, 0.0);
	xx[nn - 1] = yy[nn - 1] / uu[nn - 1][nn - 1];
	for (i = nn - 2; ii >= 0; ii--)
	{
		double sum4 = 0;
		for (kk = ii + 1; kk < nn; kk++)
			sum4 += uu[ii][kk] * xx[kk];
		xx[ii] = (yy[ii] - sum4) / uu[ii][ii];
	}
	vector<double> u_v_r(nn+1,0.0);
	u_v_r[0]=0.0;
	for(int i=1;i<nn+1;i++)
		u_v_r[i]=xx[i-1];
	//for(int i =0;i<n+1;i++)
		//cout<<u_v[i]<<" "<<endl;
	double ansur = u_v_r[ansNode1] - u_v_r[ansNode2];
	double ansr = ansur/input_i;

	
	cout << "��Ч��ѹΪ:" << ansu << endl
		<< "��Ч����Ϊ:" << ansr << endl;
	

	system("pause");

	return 0;
}
