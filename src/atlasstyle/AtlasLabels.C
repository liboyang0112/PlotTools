#include "atlasstyle/AtlasLabels.h"

#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TPad.h"
#include "TMarker.h"


void ATLASLabel(Double_t x,Double_t y,const char* text,Color_t color, const char* lumi, const char* analysis, const char* channelname) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
//   l.SetTextFont(72);
  l.SetTextFont(73);
  l.SetTextSize(15);
  l.SetTextColor(color);

  double delx = 0.06*696*gPad->GetWh()/(472*gPad->GetWw());
  double dely = 0.05;
  l.DrawLatex(x,y,"ATLAS");
  if (text) {
    TLatex p;
    p.SetNDC();
//     p.SetTextFont(42);
    p.SetTextFont(43);
    p.SetTextSize(15);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
    if(channelname){
      p.DrawLatex(x,y-dely,lumi);
      p.DrawLatex(x,y-2*dely,analysis);
      p.SetTextSize(23);
      p.DrawLatex(x,y-3*dely-0.02,channelname);
    }
  }
}

void ATLASLabelOld(Double_t x,Double_t y,bool Preliminary,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
  if (Preliminary) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+0.115,y,"Preliminary");
  }
}

void ATLASLabelNew(Double_t x,Double_t y, const char* text, Color_t color, float text_size)
{

  TLatex l;
  l.SetNDC();
  l.SetTextFont(73);
  l.SetTextColor(color);
  l.SetTextSize(text_size);

  double delx = 0.16;

  l.DrawLatex(x,y,"ATLAS");
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(43);
    p.SetTextColor(color);
    p.SetTextSize(text_size);
    p.DrawLatex(x+delx,y,text);
  }
  return;
}

void ATLASVersion(const char* version,Double_t x,Double_t y,Color_t color) 
{

  if (version) {
    char versionString[100];
    sprintf(versionString,"Version %s",version);
    TLatex l;
    l.SetTextAlign(22); 
    l.SetTextSize(0.04); 
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextColor(color);
    l.DrawLatex(x,y,versionString);
  }
}

