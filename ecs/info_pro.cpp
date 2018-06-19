#include "info_pro.h"
#include "LMatrix.h"
#include <math.h>
#include <iomanip>      
#include <ctime>        
#include <time.h> 
#include <algorithm>
const float noiseX = 4;  //噪点系数
//lwlr的参数
const float SmoothnessCoefficient = 2.55;
const float lAddX = 1.1;
//三次指数参数
const float tAddX = 1.1;
const float a = 0.585;  //指数系数
//cig_lwlr跨度
const int cig_span = 6;
const float cig_SmoothnessCoefficient = 5;
//acl_lwlr
const float acl_SmoothnessCoefficient = 0.14;
//gray_lwlr
const float gray_SmoothnessCoefficient = 0.6;


int predict_daySpan, history_daySpan; //需要预测的时间跨度和总的历史时间跨度
vector<int> vflavors_pridict_nums;  //整理好的预测值标签
string dim, predict_begin_time, predict_end_time, history_begin_time, history_end_time;
vector< vector<int> > sequences;
time_t predict_begin_time_t, predict_end_time_t, history_begin_time_t, history_end_time_t;
PhyInfo phyinfo;
FlavorsInfo flavorsinfo;

void get_info(char * info[MAX_INFO_NUM])
{
	int num_of_vm;

	string line = info[0];
	int firstSpace = line.find_first_of(' ');
	int secondSpace = line.find_last_of(' ');
	phyinfo.phycpu = atoi(line.substr(0, firstSpace).c_str());
	phyinfo.phymem = atoi(line.substr(firstSpace + 1, secondSpace - firstSpace - 1).c_str()) * 1024; 
	phyinfo.phyhard = atoi(line.substr(secondSpace + 1).c_str());

	line = info[2];
	num_of_vm = atoi(line.c_str());

	for (int i = 3; i < 3 + num_of_vm; i++)
	{
		line = info[i];
		int firstSpace = line.find_first_of(' ');
		int secondSpace = line.find_last_of(' ');
		flavorsinfo.vflavors.push_back(line.substr(0, firstSpace));
		flavorsinfo.vcpus.push_back(atoi(line.substr(firstSpace + 1, secondSpace - firstSpace - 1).c_str()));
		flavorsinfo.vmems.push_back(atoi(line.substr(secondSpace + 1).c_str()));
	}

	dim = info[3 + num_of_vm + 1];
	predict_begin_time = info[3 + num_of_vm + 3];
	predict_end_time = info[3 + num_of_vm + 4];
	predict_begin_time_t = str2time(predict_begin_time.c_str(), "%Y-%m-%d %H:%M:%S");
	predict_end_time_t = str2time(predict_end_time.c_str(), "%Y-%m-%d %H:%M:%S");
	predict_daySpan = ceil((float)(predict_end_time_t - predict_begin_time_t) / oneDayLong);

	//对flavors排序,flavor1~flavor15  加了这一步骤,对结果没影响
	for (int i = 0; i < flavorsinfo.vflavors.size() - 1; i++)
	{
		for (int j = i + 1; j < flavorsinfo.vflavors.size(); j++)
		{
			if (atoi(flavorsinfo.vflavors[j].substr(sizeof("flavor") - 1).c_str()) < atoi(flavorsinfo.vflavors[i].substr(sizeof("flavor") - 1).c_str()))
			{
				//交换
				string vflavor_tmp;
				int vcpu_tmp, vmem_temp;

				vflavor_tmp = flavorsinfo.vflavors[i];
				vcpu_tmp = flavorsinfo.vcpus[i];
				vmem_temp = flavorsinfo.vmems[i];

				flavorsinfo.vflavors[i] = flavorsinfo.vflavors[j];
				flavorsinfo.vcpus[i] = flavorsinfo.vcpus[j];
				flavorsinfo.vmems[i] = flavorsinfo.vmems[j];

				flavorsinfo.vflavors[j] = vflavor_tmp;
				flavorsinfo.vcpus[j] = vcpu_tmp;
				flavorsinfo.vmems[j] = vmem_temp;
			}
		}
	}

	cout << "phycpu: " << phyinfo.phycpu << endl;
	cout << "phymem: " << phyinfo.phymem << endl;
	cout << "phyhard: " << phyinfo.phyhard << endl;
	cout << "num_of_vm: " << num_of_vm << endl;
	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << " " << flavorsinfo.vcpus[i] << " " << flavorsinfo.vmems[i] << endl;
	}
	cout << "dim: " << dim << endl;
	cout << "predict_begin_time: " << predict_begin_time << endl;
	cout << "predict_end_time: " << predict_end_time << endl;
	cout << "predict_daySpan: " << predict_daySpan << endl;
}
void get_data(char * data[MAX_DATA_NUM], int data_num)
{
	string line = data[0];
	int firstTab = line.find_first_of('\t');
	int secondTab = line.find_last_of('\t');
	int lastSpace = line.find_last_of(' ');
	history_begin_time = line.substr(secondTab + sizeof("\t") - 1, lastSpace - secondTab - 1);
	history_begin_time.append(" 00:00:00");
	cout << history_begin_time << endl;
	line = data[data_num - 1];
	firstTab = line.find_first_of('\t');
	secondTab = line.find_last_of('\t');
	lastSpace = line.find_last_of(' ');
	history_end_time = line.substr(secondTab + sizeof("\t") - 1, lastSpace - secondTab - 1);
	history_end_time.append(" 23:59:59");
	cout << history_end_time << endl;

	history_begin_time_t = str2time(history_begin_time.c_str(), "%Y-%m-%d %H:%M:%S");
	cout << history_begin_time_t << endl;
	history_end_time_t = str2time(history_end_time.c_str(), "%Y-%m-%d %H:%M:%S");
	cout << history_end_time_t << endl;
	history_daySpan = ceil((float)(history_end_time_t - history_begin_time_t) / oneDayLong);
	cout << "history_daySpan: " << history_daySpan << endl;

	/*// get sequences
	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		vector<int> sequence(history_daySpan, 0); //history_daySpan ints with value 0;
												  //search data
		for (int j = 0; j < data_num; j++)
		{
			string line = data[j];
			int firstTab = line.find_first_of('\t');
			int secondTab = line.find_last_of('\t');
			int lastSpace = line.find_last_of(' ');
			string vflavor = line.substr(firstTab + sizeof("\t") - 1, secondTab - firstTab - sizeof("\t") + 1);
			time_t time = str2time(line.substr(secondTab + sizeof("\t") - 1, lastSpace - secondTab - 1).append(" 00:00:00").c_str(), "%Y-%m-%d %H:%M:%S");
			int dayDiff = (time - history_begin_time_t) / oneDayLong; // used as sequence index
																	  //            cout << vflavor << "\t" << time << "\t" << dayDiff << endl;
			if (vflavor == flavorsinfo.vflavors[i]) //this kind of flavor need to be predicted
			{
				sequence[dayDiff]++;   //add a record at daydiff index
			}
		}
		sequences.push_back(sequence);
	}*/
	vector<int> sequence(history_daySpan, 0);

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		sequences.push_back(sequence);
	}

	for (int j = 0; j < data_num; j++)
	{

		string line = data[j];
		int firstTab = line.find_first_of('\t');
		int secondTab = line.find_last_of('\t');
		int lastSpace = line.find_last_of(' ');
		string vflavor = line.substr(firstTab + sizeof("\t") - 1, secondTab - firstTab - sizeof("\t") + 1);
		time_t time = str2time(line.substr(secondTab + sizeof("\t") - 1, lastSpace - secondTab - 1).append(" 00:00:00").c_str(), "%Y-%m-%d %H:%M:%S");
		int dayDiff = (time - history_begin_time_t) / oneDayLong; // used as sequence index
		vector<string>::iterator vf_post = find(flavorsinfo.vflavors.begin(), flavorsinfo.vflavors.end(), vflavor);
		if (vf_post != flavorsinfo.vflavors.end())
		{
			int dis = distance(flavorsinfo.vflavors.begin(), vf_post);//cout << vflavor << "\t" << time << "\t" << dayDiff << endl;
			sequences[dis][dayDiff]++;
		}

	}
	//去燥
	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		float aveValue = accumulate(sequences[i].begin(),sequences[i].end(),0) / history_daySpan+ 1;
		for (int k = 0; k<history_daySpan; k++)
		{
			if (sequences[i][k]>noiseX * aveValue)
				sequences[i][k] = noiseX * aveValue;
		}
	}
	//

	for (int i = 0; i < sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		int sum = 0;
		for (int j = 0; j < sequences[i].size(); j++)
		{
			sum += sequences[i][j];
			cout << j << ": " << sequences[i][j] << ", ";
		}
		cout << "total_num: " << sum << endl;
		cout << endl;
	}
}

time_t str2time(string str, string format)
{
	struct tm *tmp_time = (struct tm*)malloc(sizeof(struct tm));
	strptime(str.c_str(), format.c_str(), tmp_time);
	time_t t = mktime(tmp_time);
	free(tmp_time);
	return t;

	/*tm mytm = {};
	istringstream iss(str);
	iss >> std::get_time(&mytm, format.c_str());
	time_t t = mktime(&mytm);
	return t;*/
}


void simple_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)

		for (int j = 0; j < sequences[i].size() - predict_daySpan + 1; j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + j, sequences[i].begin() + j + predict_daySpan, 0));
		}
		daySpan_sequences.push_back(daySpan_sequence);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
	}

	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		vector<int> temp = daySpan_sequences[i];  //temp for daySpan_sequences
		for (int j = 0; j < predict_daySpan; j++)
		{
			int sum = accumulate(temp.begin(), temp.end(), 0);
			float average = (float)sum / temp.size();
			float last_daySpan = temp.back();


			float w_last_daySpan = 0.95;
			float w_average = 1 - w_last_daySpan;
			int predict_result = round(w_average * average + w_last_daySpan * last_daySpan);
			//            cout << predict_result;
			temp.push_back(predict_result);
		}
		//vflavors_pridict_nums.push_back(temp.back()); //再次整理数据集
		temp.back()<0 ? vflavors_pridict_nums.push_back(0) : vflavors_pridict_nums.push_back(temp.back());
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}

void lwlr_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector<LMatrix<float>> Y;//各类型虚拟机的销售标签数据
	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		LMatrix<float>y(int(sequences[i].size() / predict_daySpan), 1);

		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0));
			y[j][0] = accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0);
		}
		daySpan_sequences.push_back(daySpan_sequence);
		Y.push_back(y);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << Y[i][j][0] << ", ";
		}
		cout << endl;
	}

	LMatrix<float> x(Y[0].RowLen, 2, 1);  //特征矩阵初始化
	for (int i = 0; i < x.RowLen; i++)
		x[i][1] = i;
	LMatrix<float> testPoint(1, 2, 1);
	testPoint[0][1] = x.RowLen;
	LMatrix<float> diffmat, tempmat;

	for (int i = 0; i < Y.size(); i++)
	{
		LMatrix<float> weight;     //权重矩阵
		weight.eye(x.RowLen);
		for (int i = 0; i < x.RowLen; i++)//计算权重
		{
			diffmat = testPoint - x.GetRow(i);
			tempmat = diffmat * diffmat.T();
			float tem = tempmat[0][0] * 1;
			float we = exp(tem / (-2.0*SmoothnessCoefficient*SmoothnessCoefficient));
			weight[i][i] = we;
		}
		LMatrix<float> xTx, ws;
		xTx = x.T()*(weight*x);
		ws = xTx.INV() * (x.T() * (weight * Y[i]));//计算出回归系数的一个估计
		cout << (testPoint * ws)[0][0] << endl;
		//vflavors_pridict_nums.push_back((testPoint * ws)[0][0]); //再次整理数据集
		(testPoint * ws)[0][0]<0 ? vflavors_pridict_nums.push_back(1) : vflavors_pridict_nums.push_back((testPoint * ws)[0][0]*lAddX);
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}

void cig_lwlr_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector<LMatrix<float>> Y;//各类型虚拟机的销售标签数据
	vector<LMatrix<float>> X;
	vector<LMatrix<float>> TESTPOINT;
	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		LMatrix<float>y(int(sequences[i].size() / predict_daySpan) - cig_span + 1, 1);
		LMatrix<float>x(int(sequences[i].size() / predict_daySpan) - cig_span + 1, cig_span - 1);
		LMatrix<float> testPoint(1, cig_span - 1);

		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0));
		}
		for (int j = 0; j < int(sequences[i].size() / predict_daySpan) - cig_span + 1; j++)
		{
			y[j][0] = float(daySpan_sequence[j + cig_span - 1]);
			for (int k = 0; k < cig_span - 1; k++)
			{
				x[j][k] = daySpan_sequence[j + k];
			}
		}
		for (int j = 0; j < cig_span - 1; j++)
		{
			testPoint[0][j] = daySpan_sequence[int(sequences[i].size() / predict_daySpan) - cig_span + j];
		}
		daySpan_sequences.push_back(daySpan_sequence);
		Y.push_back(y);
		X.push_back(x);
		TESTPOINT.push_back(testPoint);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
	}


	LMatrix<float> diffmat, tempmat;

	for (int i = 0; i < Y.size(); i++)
	{
		LMatrix<float> weight;     //权重矩阵
		weight.eye(X[i].RowLen);
		for (int j = 0; j < X[i].RowLen; j++)//计算权重
		{
			diffmat = TESTPOINT[i] - X[i].GetRow(j);
			tempmat = diffmat * diffmat.T();
			float tem = tempmat[0][0] * 1;
			float we = exp(tem / (-2.0*cig_SmoothnessCoefficient*cig_SmoothnessCoefficient));

			weight[j][j] = we;
			//weight.show();
		}
		LMatrix<float> xTx, ws;
		xTx = X[i].T()*(weight*X[i]);
		ws = xTx.INV() * (X[i].T() * (weight * Y[i]));//计算出回归系数的一个估计
		cout << (TESTPOINT[i] * ws)[0][0] << endl;
		//vflavors_pridict_nums.push_back((TESTPOINT[i] * ws)[0][0]); //再次整理数据集
		(TESTPOINT[i] * ws)[0][0] < 0 ? vflavors_pridict_nums.push_back(0) : vflavors_pridict_nums.push_back((TESTPOINT[i] * ws)[0][0]);
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}

void aclwlr_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector<LMatrix<float>> Y;//各类型虚拟机的销售标签数据
	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		LMatrix<float>y(sequences[0].size() - predict_daySpan + 1, 1);

		for (int j = 0; j < sequences[i].size() - predict_daySpan + 1; j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + j, sequences[i].begin() + j + predict_daySpan, 0));
			y[j][0] = accumulate(sequences[i].begin() + j, sequences[i].begin() + j + predict_daySpan, 0);
		}
		daySpan_sequences.push_back(daySpan_sequence);
		Y.push_back(y);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
	}

	LMatrix<float> x(Y[0].RowLen, 2, 1);  //特征矩阵初始化
	for (int i = 0; i < x.RowLen; i++)
		x[i][1] = i;
	LMatrix<float> testPoint(1, 2, 1);
	testPoint[0][1] = x.RowLen;
	LMatrix<float> diffmat, tempmat;

	for (int i = 0; i < Y.size(); i++)
	{
		LMatrix<float> weight;     //权重矩阵
		weight.eye(x.RowLen);
		for (int i = 0; i < x.RowLen; i++)//计算权重
		{
			diffmat = testPoint - x.GetRow(i);
			tempmat = diffmat * diffmat.T();
			float tem = tempmat[0][0] * 0.001;
			float we = exp(tem / (-2.0*acl_SmoothnessCoefficient*acl_SmoothnessCoefficient));
			weight[i][i] = we;
		}
		LMatrix<float> xTx, ws;
		xTx = x.T()*(weight*x);
		ws = xTx.INV() * (x.T() * (weight * Y[i]));//计算出回归系数的一个估计
		//vflavors_pridict_nums.push_back((testPoint * ws)[0][0]); //再次整理数据集
		(testPoint * ws)[0][0]<0 ? vflavors_pridict_nums.push_back(1) : vflavors_pridict_nums.push_back((testPoint * ws)[0][0]);
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}

void gray_lwlr_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector<LMatrix<float>> Y;//各类型虚拟机的销售标签数据
	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		LMatrix<float>y(int(sequences[i].size() / predict_daySpan), 1);

		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0));
			y[j][0] = accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0);
		}
		daySpan_sequences.push_back(daySpan_sequence);
		Y.push_back(y);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << Y[i][j][0] << ", ";
		}
		cout << endl;
	}

	LMatrix<float> x(Y[0].RowLen, 2, 1);  //特征矩阵初始化
	for (int i = 0; i < x.RowLen; i++)
		x[i][1] = i;
	LMatrix<float> testPoint(1, 2, 1);
	testPoint[0][1] = x.RowLen;
	LMatrix<float> diffmat, tempmat;

	for (int i = 0; i < Y.size(); i++)
	{
		LMatrix<float> weight;     //权重矩阵
		weight.eye(x.RowLen);
		for (int i = 0; i < x.RowLen; i++)//计算权重
		{
			diffmat = testPoint - x.GetRow(i);
			tempmat = diffmat * diffmat.T();
			float tem = tempmat[0][0] * 1;
			float we = exp(tem / (-2.0*gray_SmoothnessCoefficient*gray_SmoothnessCoefficient));
			weight[i][i] = we;
		}
		LMatrix<float> xTx, ws;
		xTx = x.T()*(weight*x);
		ws = xTx.INV() * (x.T() * (weight * Y[i]));//计算出回归系数的一个估计
		cout << (testPoint * ws)[0][0] << endl;
		//vflavors_pridict_nums.push_back((testPoint * ws)[0][0] - Y[i][x.RowLen - 1][0]); //再次整理数据集
		(testPoint * ws)[0][0] - Y[i][x.RowLen - 1][0] < 0 ? vflavors_pridict_nums.push_back(0) : vflavors_pridict_nums.push_back((testPoint * ws)[0][0] - Y[i][x.RowLen - 1][0]);
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}






void ThreeIndex_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector< vector<float> > S0;//各类型虚拟机的销售数据
	vector< vector<float> > S1;
	vector< vector<float> > S2;
	vector< vector<float> > S3;

	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		vector<float>s0;
		vector<float>initdata;
		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0));
			s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
		}
		initdata.push_back((s0[0] + s0[1] + s0[2]) / 3);
		daySpan_sequences.push_back(daySpan_sequence);
		S0.push_back(s0);
		S1.push_back(initdata);
		S2.push_back(initdata);
		S3.push_back(initdata);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << S0[i][j] << ", ";
		}
		cout << endl;
	}
	//计算S1
	for (int i = 0; i < S0.size(); i++)
	{
		for (int j = 1; j < S0[i].size(); j++)
		{
			S1[i].push_back(a*S0[i][j] + (1 - a)*S1[i][j - 1]);
		}
	}
	//计算S2
	for (int i = 0; i < S1.size(); i++)
	{
		for (int j = 1; j < S1[i].size(); j++)
		{
			S2[i].push_back(a*S1[i][j] + (1 - a)*S2[i][j - 1]);
		}
	}
	//计算S3
	for (int i = 0; i < S2.size(); i++)
	{
		for (int j = 1; j < S2[i].size(); j++)
		{
			S3[i].push_back(a*S2[i][j] + (1 - a)*S3[i][j - 1]);
		}
		cout << S3[i].back() << endl;
	}

	//开始预测
	for (int i = 0; i < S0.size(); i++)
	{
		float At, Bt, Ct;
		At = 3 * S1[i].back() - 3 * S2[i].back() + S3[i].back();
		Bt = (a / (2 * pow(1 - a, 2)))*((6 - 5 * a)*S1[i].back() - 2 * (5 - 4 * a)*S2[i].back() + (4 - 3 * a)*S3[i].back());
		Ct = (a*a / (2 * pow(1 - a, 2)))*(S1[i].back() - 2 * S2[i].back() + S3[i].back());
		cout << At + Bt + Ct << endl;
		(At + Bt + Ct - S0[i].back()) <= 0 ? vflavors_pridict_nums.push_back(1) : vflavors_pridict_nums.push_back((At + Bt + Ct - S0[i].back())*tAddX);
		//vflavors_pridict_nums.push_back(At+Bt+Ct-S0[i].back()); //预测下一期
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}















/*
void cubicSmooth5(vector<int > in, vector< int > &out)
{
	
	
		
		int N = in.size();
		
			int i;
			if (N < 5)
			{
				
				for (i = 0; i <= N - 1; i++)
					out.push_back(in[i]);
			}

			else
			{
				
				out.push_back((69.0 * in[0] + 4.0 * in[1] - 6.0 * in[2] + 4.0 * in[3] - in[4]) / 70.0);
				
				out.push_back((2.0 * in[0] + 27.0 * in[1] + 12.0 * in[2] - 8.0 * in[3] + 2.0 * in[4]) / 35.0);
				for (i = 2; i <= N - 3; i++)
				{
					out.push_back((-3.0 * (in[i - 2] + in[i + 2]) + 12.0 * (in[i - 1] + in[i + 1]) + 17.0 * in[i]) / 35.0);
				}
				out.push_back((2.0 * in[N - 5] - 8.0 * in[N - 4] + 12.0 * in[N - 3] + 27.0 * in[N - 2] + 2.0 * in[N - 1]) / 35.0);
				out.push_back((-in[N - 5] + 4.0 * in[N - 4] - 6.0 * in[N - 3] + 4.0 * in[N - 2] + 69.0 * in[N - 1]) / 70.0);
			}
			
			
	
	return;
}
void cubicSmooth7(vector<vector<int >> in, vector<vector< int >> &out)
{
	for (int q = 0; q < in.size(); q++)
	{
		vector<int > seq;
		int N = in[q].size();
		int i;
		if (N < 7)
		{
			for (i = 0; i <= N - 1; i++)
			{
				seq.push_back(in[q][i]);
			}
		}
		else
		{
			seq.push_back((39.0 * in[q][0] + 8.0 * in[q][1] - 4.0 * in[q][2] - 4.0 * in[q][3] +
				1.0 * in[q][4] + 4.0 * in[q][5] - 2.0 * in[q][6]) / 42.0);
			seq.push_back((8.0 * in[q][0] + 19.0 * in[q][1] + 16.0 * in[q][2] + 6.0 * in[q][3] -
				4.0 * in[q][4] - 7.0* in[q][5] + 4.0 * in[q][6]) / 42.0);
			seq.push_back((-4.0 * in[q][0] + 16.0 * in[q][1] + 19.0 * in[q][2] + 12.0 * in[q][3] +
				2.0 * in[q][4] - 4.0 * in[q][5] + 1.0 * in[q][6]) / 42.0);
			for (i = 3; i <= N - 4; i++)
			{
				seq.push_back((-2.0 * (in[q][i - 3] + in[q][i + 3]) +
					3.0 * (in[q][i - 2] + in[q][i + 2]) +
					6.0 * (in[q][i - 1] + in[q][i + 1]) + 7.0 * in[q][i]) / 21.0);
			}
			seq.push_back((-4.0 * in[q][N - 1] + 16.0 * in[q][N - 2] + 19.0 * in[q][N - 3] +
				12.0 * in[q][N - 4] + 2.0 * in[q][N - 5] - 4.0 * in[q][N - 6] + 1.0 * in[q][N - 7]) / 42.0);
			seq.push_back((8.0 * in[q][N - 1] + 19.0 * in[q][N - 2] + 16.0 * in[q][N - 3] +
				6.0 * in[q][N - 4] - 4.0 * in[q][N - 5] - 7.0 * in[q][N - 6] + 4.0 * in[q][N - 7]) / 42.0);
			seq.push_back((39.0 * in[q][N - 1] + 8.0 * in[q][N - 2] - 4.0 * in[q][N - 3] -
				4.0 * in[q][N - 4] + 1.0 * in[q][N - 5] + 4.0 * in[q][N - 6] - 2.0 * in[q][N - 7]) / 42.0);
		}
		out.push_back(seq);
	}
}
*/
/*
void ThreeIndex_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector< vector<float> > S0;//各类型虚拟机的销售数据
	vector< vector<float> > S1;
	vector< vector<float> > S2;
	vector< vector<float> > S3;

	
	
	for (int i = 0; i < sequences.size(); i++)
	{
		
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		vector<float>s0;
		vector<float>initdata;

		vector<vector<int >>sequences_copy = sequences;
		sequences.clear();
		cubicSmooth7(sequences_copy,sequences);
		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0));
			s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
		}
	
		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.begin()+j, 0));
		}
		initdata.push_back((s0[0] + s0[1] + s0[2]) / 3);
		daySpan_sequences.push_back(daySpan_sequence);
		S0.push_back(s0);
		S1.push_back(initdata);
		S2.push_back(initdata);
		S3.push_back(initdata);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << S0[i][j] << ", ";
		}
		cout << endl;
	}
	//计算S1
	for (int i = 0; i < S0.size(); i++)
	{
		for (int j = 1; j < S0[i].size(); j++)
		{
			S1[i].push_back(a*S0[i][j] + (1 - a)*S1[i][j - 1]);
		}
	}
	//计算S2
	for (int i = 0; i < S1.size(); i++)
	{
		for (int j = 1; j < S1[i].size(); j++)
		{
			S2[i].push_back(a*S1[i][j] + (1 - a)*S2[i][j - 1]);
		}
	}
	//计算S3
	for (int i = 0; i < S2.size(); i++)
	{
		for (int j = 1; j < S2[i].size(); j++)
		{
			S3[i].push_back(a*S2[i][j] + (1 - a)*S3[i][j - 1]);
		}
		cout << S3[i].back() << endl;
	}
	//开始预测
	for (int i = 0; i < S0.size(); i++)
	{
		float At, Bt, Ct;
		At = 3 * S1[i].back() - 3 * S2[i].back() + S3[i].back();
		Bt = (a / (2 * pow(1 - a, 2)))*((6 - 5 * a)*S1[i].back() - 2 * (5 - 4 * a)*S2[i].back() + (4 - 3 * a)*S3[i].back());
		Ct = (a*a / (2 * pow(1 - a, 2)))*(S1[i].back() - 2 * S2[i].back() + S3[i].back());
		cout << At + Bt + Ct << endl;
		At + Bt + Ct - S0[i].back() < 0 ? vflavors_pridict_nums.push_back(0) : vflavors_pridict_nums.push_back(At + Bt + Ct - S0[i].back());

		//vflavors_pridict_nums.push_back(At + Bt + Ct - S0[i].back()); //预测下一期
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i] << endl;
	}
}
*/


/*
void ThreeIndex_predict()
{
	// init sequences
	vector< vector<int> > daySpan_sequences;
	vector< vector<float> > S0;//各类型虚拟机的销售数据
	vector< vector<float> > S1;
	vector< vector<float> > S2;
	vector< vector<float> > S3;
	for (int i = 0; i < sequences.size(); i++)
	{
		vector<int> daySpan_sequence;  //按时间跨度整理好的数据集(初始数据)
		vector<float>s0;
		vector<float>initdata;
		bool issuit = false;
		for (int j = 0; j < int(sequences[i].size() / predict_daySpan); j++)
		{
			
			if ((j + 1)< int(sequences[i].size() / predict_daySpan) && abs(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
				sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0)- 
				accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan,
					sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 2) * predict_daySpan, 0))>=15)
			{
				if (accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
					sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0) >
					accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan,
						sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 2) * predict_daySpan, 0))
				{
					daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
						sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0)-2
					);
					s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
					daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j+1) * predict_daySpan,
						sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j +2) * predict_daySpan, 0) + 0);
					s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
					issuit = true;
				}
			else 
			{
				daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
					sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0) +0
				);
				s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
				daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j+1) * predict_daySpan,
					sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j +2) * predict_daySpan, 0) -2);
				s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
				issuit = true;
			}
			}
			else
			{
				if (issuit)
				{
				
					issuit = false;
					
				}
				else 
				{
					
					daySpan_sequence.push_back(accumulate(sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + j * predict_daySpan,
					sequences[i].begin() + int(sequences[i].size() % predict_daySpan) + (j + 1) * predict_daySpan, 0));
					s0.push_back(accumulate(daySpan_sequence.begin(), daySpan_sequence.end(), 0));
				}
			}
			
		}
		initdata.push_back((s0[0] + s0[1] + s0[2]) / 3);
		daySpan_sequences.push_back(daySpan_sequence);
		S0.push_back(s0);
		S1.push_back(initdata);
		S2.push_back(initdata);
		S3.push_back(initdata);
	}
	cout << "daySpan_sequences: " << endl;
	for (int i = 0; i < daySpan_sequences.size(); i++)
	{
		cout << flavorsinfo.vflavors[i] << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << daySpan_sequences[i][j] << ", ";
		}
		cout << endl;
		for (int j = 0; j < daySpan_sequences[i].size(); j++)
		{
			cout << i * daySpan_sequences.size() + j << ":" << S0[i][j] << ", ";
		}
		cout << endl;
	}
	//计算S1
	for (int i = 0; i < S0.size(); i++)
	{
		for (int j = 1; j < S0[i].size(); j++)
		{
			S1[i].push_back(a*S0[i][j] + (1 - a)*S1[i][j - 1]);
		}
	}
	//计算S2
	for (int i = 0; i < S1.size(); i++)
	{
		for (int j = 1; j < S1[i].size(); j++)
		{
			S2[i].push_back(a*S1[i][j] + (1 - a)*S2[i][j - 1]);
		}
	}
	//计算S3
	for (int i = 0; i < S2.size(); i++)
	{
		for (int j = 1; j < S2[i].size(); j++)
		{
			S3[i].push_back(a*S2[i][j] + (1 - a)*S3[i][j - 1]);
		}
		cout << S3[i].back() << endl;
	}
	//开始预测

	for (int i = 0; i < S0.size(); i++)
	{
		float At, Bt, Ct;
		At = 3 * S1[i].back() - 3 * S2[i].back() + S3[i].back();
		Bt = (a / (2 * pow(1 - a, 2)))*((6 - 5 * a)*S1[i].back() - 2 * (5 - 4 * a)*S2[i].back() + (4 - 3 * a)*S3[i].back());
		Ct = (a*a / (2 * pow(1 - a, 2)))*(S1[i].back() - 2 * S2[i].back() + S3[i].back());
		cout << At + Bt + Ct << endl;
		At + Bt + Ct - S0[i].back() < 0 ? vflavors_pridict_nums.push_back(0) : vflavors_pridict_nums.push_back(At + Bt + Ct - S0[i].back());
	
		//vflavors_pridict_nums.push_back(At + Bt + Ct - S0[i].back()); //预测下一期
	}

	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
	
		cout << flavorsinfo.vflavors[i] << ":\t" << vflavors_pridict_nums[i]<< endl;
	}
}*/
