#include "LatexChart.h"
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
using namespace std;

void LatexChart::set(std::string row, std::string column, float nominal, float errorup, float errordown){
	observable obs(nominal,errorup,errordown);
	set(row,column,obs);
}

void LatexChart::set(std::string row, std::string column, double nominal, double errorup, double errordown){
	observable obs(nominal,errorup,errordown);
	set(row,column,obs);
}

void LatexChart::set(std::string row, std::string column, observable obs){
	auto iter = find(rows.begin(),rows.end(),row);
	if (iter == rows.end()) rows.push_back(row);
	iter = find(columns.begin(),columns.end(),column);
	if (iter == columns.end()) columns.push_back(column);
	content[row][column].push_back(obs);
}

void LatexChart::clear(){
	rows.clear();
	columns.clear();
	content.clear();
}

void LatexChart::reset(){
	for(auto row : content)
		for(auto column : row.second)
			column.second.clear();
}

void LatexChart::print(std::string filename){
	printf("LatexChart::print() : Print to file: %s\n",filename.c_str());
	ofstream *file = new ofstream();
	(*file).open(filename+".tex");
	(*file)<<"\\centering\n";
	int ncolumn = columns.size();
	int currentnrow = 0;
	if(ncolumn > maxcolumn){
		int nvec = ncolumn/maxcolumn;
		if(ncolumn%maxcolumn) nvec+=1;
		int nlong = ncolumn%nvec;
		int averagelow = ncolumn/nvec;
		if(!ncolumn%maxcolumn) averagelow-=1;
		vector<string> new_columns;
		int count = 0;
		char nfile = '0';
		for (int ivec = 0; ivec < nvec; ++ivec)
		{
			if((currentnrow+=rows.size()+1) > maxrow){
				currentnrow=rows.size()+1;
				file->close();
				(*file).open(filename+"_"+char(++nfile)+".tex");
				(*file)<<"\\centering\n";
			}
			for (int i = 0; i < averagelow+1; ++i)
			{
				if(i == averagelow && ivec >= nlong) continue;
				new_columns.push_back(columns[count]);
				count ++;
			}
			writeContent(new_columns, file);
			new_columns.clear();
		}
	}else{
		writeContent(columns, file);
	}
	file->close();
}

void LatexChart::writeContent(std::vector<std::string> new_columns, std::ofstream* file){
	(*file)<<"\\begin{tabular}{|";
	for(auto new_column: new_columns) (*file)<<"c|";
	if(!threelinetable) (*file)<<"c|} \\hline\n";
	else (*file)<<"c|} \\toprule\\toprule\n";
	//==============================column title=====================================
	for(auto new_column: new_columns) (*file)<<" & "<<new_column;
	if(!threelinetable) (*file)<<"\\\\\\hline\n";
	else "\\\\\\midrule\n";
	//==============================table content=====================================
	(*file)<<fixed<<setprecision(2);
	for(auto row: rows){
		(*file)<<row;
		for(auto new_column: new_columns) {
			(*file)<<" & ";
			if(content[row].find(new_column) == content[row].end())
				(*file)<<" /";
			else{
				for (int icontent = 0; icontent < content[row][new_column].size(); ++icontent)
				{
					observable* target = grabContent(row,new_column,icontent);
					if(!target) {
						break;
					}
					(*file)<<"$"<<target->nominal;
					if(target->error) {
						if(target->error == target->errordown)
							(*file)<<"\\pm"<<target->error;
						else
							(*file)<<"^{+"<<target->error<<"}_{-"<<target->errordown<<"}";
					}
					(*file)<<"$";
					if(icontent != content[row][new_column].size()-1) (*file)<<" / ";
				}
			}
		}
		if(!threelinetable) (*file)<<"\\\\\\hline\n";
		else (*file)<<"\\\\\n";
	}
	if(threelinetable) (*file)<<"\\bottomrule\\bottomrule\\\\\n";
	(*file)<<"\\end{tabular}\n";
}

observable* LatexChart::grabContent(string _row, string _colume, int _icontent){
  for(auto &row : content){
    if(row.first != _row) continue;
    for(auto &column : row.second){
      if(column.first != _colume) continue;
      if(_icontent >= column.second.size()) {
        if(debug) printf("LatexChart::grabContent()  WARNING: table content doesn't exist: row = %s. column = %s, icontent = %d\n", _row.c_str(), _colume.c_str(), _icontent);
      }
      else return &(column.second.at(_icontent));
      if(debug) printf("LatexChart::grabContent()  WARNING: table column doesn't exist: row = %s. column = %s\n", _row.c_str(), _colume.c_str());
      return 0;
    }
    if(debug) printf("LatexChart::grabContent()  WARNING: table row doesn't exist: row = %s.\n", _row.c_str());
    return 0;
  }
  return 0;
}

void LatexChart::add(LatexChart *target){
	for(auto row: rows){
		for(auto column: columns){
			for (int icontent = 0; icontent < content[row][column].size(); icontent++)
			{
				observable* thiscontent = grabContent(row,column,icontent);
				observable* targetcontent = target->grabContent(row,column,icontent);
				if(thiscontent && targetcontent) *thiscontent += *targetcontent;
			}
		}
	}
}

void LatexChart::concate(LatexChart *target){
	for(auto row: rows){
		for(auto column: columns){
			for (int icontent = 0; icontent < target->content[row][column].size(); icontent++)
			{
				observable* targetcontent = target->grabContent(row,column,icontent);
				if(targetcontent) content[row][column].push_back(*targetcontent);
			}
		}
	}
}

LatexChart* LatexChart::clone(){
	LatexChart *chart = new LatexChart(label);
	chart->maxcolumn = maxcolumn;
	chart->label = label;
	chart->caption = caption;
	chart->rows = rows;
	chart->columns = columns;
	chart->content = content;
	return chart;
};
