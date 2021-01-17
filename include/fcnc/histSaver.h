#ifndef histSaver_h
#define histSaver_h
#include <iostream>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "observable.h"

struct variable{

  TString name;
  TString title;
  int nbins;
  float xlow;
  float xhigh;
  TString unit;
  float scale;
  int rebin;
  std::vector<double>* xbins;

  variable(TString _name, TString _title, int _nbins, float _xlow, float _xhigh, TString _unit = "", float _scale = 1, int _rebin = 1,   std::vector<double>* _xbins = 0)
  :name(_name), title(_title), nbins(_nbins), xlow(_xlow), xhigh(_xhigh), unit(_unit), scale(_scale), rebin(_rebin), xbins(_xbins){}

};

struct fcncSample
{
public:
  fcncSample(TString _name, TString _title, enum EColor _color, double _norm = 1): name(_name), title(_title), color(_color), norm(_norm) {};
  TString name;
  TString title;
  double norm;
  enum EColor color;
};

class histSaver{
public:
  TString inputfilename;
  bool checkread;
  TString checkread_sample;
  TString checkread_region;
  TString checkread_variation;
  bool checkread_variable;
  int checkread_ibin;

  int nregion;
  float blinding;
  Float_t* fweight;
  Double_t* dweight;
  int weight_type;
  std::vector<TString> overlaysamples;
  TFile* inputfile;
  std::map<TString,TFile*> outputfile;
  std::map<TString,std::string> regioninTables;
  TString lumi;
  TString analysis;
  TString workflow;
  bool doROC;
  TString nominalfilename;
  TString createdNP;
  int debug;
  std::vector<variable*> v;
  Int_t nvar;
  std::vector<Float_t*> address1;
  std::vector<Double_t*> address3;
  std::vector<Int_t*> address2;
  bool dataref;
  TString trexdir;
  std::vector<TString> stackorder;
  std::map<TString,std::vector<TString>> merge_suffix;
  TString outputfilename;
  TString sensitivevariable;
  std::map<TString, std::map<TString, std::map<TString, std::vector<TH1D*> > > > plot_lib; //plot_lib[sample][region][variation][var]
  std::vector<TString> regions;
  std::vector<fcncSample> samples;
  std::vector<TString> mutedregions;
  static TFile *bufferfile;
  histSaver(TString outputfilename);
  virtual ~histSaver();
  void clearhist();
  template<typename D>
  void add(variable* v_, D* var_ = 0){
    v.push_back(v_);
    TString Dname = typeid(*var_).name();
    if(debug) printf("fill type: %s\n", Dname.Data());
    address1.push_back(0);
    address2.push_back(0);
    address3.push_back(0);
    if(var_){
      if (Dname.Contains("i")) 
        address2[v.size()-1] = (int*)var_;
      else if(Dname.Contains("f")) 
        address1[v.size()-1] = (float*)var_;
      else if(Dname.Contains("d"))
        address3[v.size()-1] = (double*)var_;
      else
        printf("unknown var type: %s\n", v_->name.Data());
    }
  }
  void add(variable* v_){v.push_back(v_);}
  void show();
  bool find_sample(TString sample);
  TH1D* grabbkghist(TString region, int ivar);
  TH1D* grabsighist(TString region, int ivar, TString signal="");
  TH1D* grabdatahist(TString region, int ivar);
  bool add_variation(TString sample,TString region,TString variation);
  void printyield(TString region);
  observable calculateYield(TString region, std::string formula, TString variation);
  double gethisterror(TH1* hist);
  observable templatesample(TString fromregion, TString variation,std::string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color,bool scaletogap, double SF = 1);
  observable templatesample(TString fromregion, TString variation,std::string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color,bool scaletogap, observable SF);
  std::vector<observable> scale_to_data(TString scaleregion, std::string formula, TString scaleVariable, std::vector<double> slices = {}, TString variation = "NOMINAL");
  void scale_sample(TString scaleregion, std::string formula, TString scaleVariable, std::vector<observable> scalefactor, std::vector<double> slices = {}, TString variation = "NOMINAL");
  int findvar(TString varname);
  std::vector<int> resolveslices(TH1D* target, const std::vector<double>* slices);
  std::map<TString,std::vector<observable>>* fit_scale_factor(std::vector<TString> *fit_regions, TString *variable, std::vector<TString> *scalesamples, const std::vector<double> *slices, TString *variation, std::vector<TString> *postfit_regions);
  std::map<TString,std::vector<observable>>* fit_scale_factor(std::vector<TString> *fit_regions, TString *variable, std::map<TString,std::map<TString,std::vector<TString>>> *scalesamples, const std::vector<double> *slices, TString *variation = 0, std::map<TString,std::map<TString,std::vector<TString>>> *postfit_regions = 0);
  void muteregion(TString region);
  void unmuteregion(TString region);
  void SetLumiAnaWorkflow(TString _lumi, TString _analysis, TString _workflow);
  void write_trexinput(TString NPname = "NOMINAL", TString writename = "", TString writeoption = "update");
  void overlay(TString _overlaysample);
  TH1D* grabhist_int(TString sample, TString region, int ivar, bool vital = 0);
  TH1D* grabhist(TString sample, TString region, TString variation, int ivar, bool vital = 0);
  TH1D* grabhist(TString sample, TString region, TString variation, TString varname, bool vital = 0);
  TH1D* grabhist(TString sample, TString region, TString varname, bool vital = 0);
  void merge_regions(TString inputregion1, TString inputregion2, TString outputregion);
  void merge_regions(std::vector<TString> inputregions, TString outputregion);
  Float_t getVal(Int_t i);
  float binwidth(int i);
  void read_sample(TString samplename, TString savehistname, TString NPname, TString sampleTitle, enum EColor color, double norm, TFile *_inputfile=0);
  void plot_stack(TString NPname,TString outputdir = ".",TString outputchartdir = ".");
  void fill_hist(TString sample, TString region, TString variation);
  void fill_hist(TString sample, TString region);
  void add_region(TString region);
  void add_sample(TString samplename, TString sampleTitle, enum EColor color);
  void init_hist(std::map<TString,std::map<TString,std::map<TString,std::vector<TH1D*>>>>::iterator sample_lib, TString region, TString variation);

  void set_weight(Float_t* _weight){ fweight = _weight; weight_type = 1;}
  void set_weight(Double_t* _weight){ dweight = _weight; weight_type = 2;}
  void write();
  // hadhad FF
  void FakeFactorMethod(TString final_region, TString _1m1lregion,TString _1l1mregion,TString _1l1nregion,TString _1n1lregion,TString _2nregion,TString variation,TString newsamplename,TString newsampletitle,std::vector<TString> tmp_regions,enum EColor color,bool SBplot);
  void FakeFactorMethod(TString final_region, TString _1m1lnmregion,TString _1lnm1mregion,TString _2nregion,TString variation,TString newsamplename,TString newsampletitle,std::vector<TString> tmp_regions,enum EColor color,bool SBplot);
  void FakeFactorMethod(TString final_region, TString _reg1mtau1ltau1b2jos,TString variation,TString newsamplename,TString newsampletitle,std::vector<TString> tmp_regions,enum EColor color,bool SBplot);
  std::vector<observable> derive_ff(std::string formula,std::vector<TString> region_numerator,std::vector<TString> region_denominator,TString variable);
  
};


#endif
