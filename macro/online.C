// display for merged GET files. Event argument is the raw GET event ID.
// Mapping deliberately matches corea/data/padHelper.hh (channels 0..63).
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <Buttons.h>
#include <TControlBar.h>
#include <TH2Poly.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TString.h>
#include <array>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

void corea_hover();

namespace CoreaOnline {
struct Frame { std::streamoff pos; uint32_t size, header, items; unsigned type, item, asad; };
std::map<uint32_t,std::vector<Frame> > eventFrames;
std::string path;
std::streamoff scanned=12;
uint32_t current=0;
int run=24;
double threshold=50;
TCanvas *canvas=nullptr;
TControlBar *controls=nullptr;
TH2Poly *charge=nullptr,*time=nullptr;
TH1D *pulse=nullptr;
std::vector<double> samples(1024*512,0);
std::vector<unsigned char> present(1024*512,0);
unsigned bestPad=1;
uint32_t be(const unsigned char *p,unsigned n) { uint32_t v=0; while(n--) v=(v<<8)|*p++; return v; }
int pad(unsigned asad,unsigned aget,unsigned ch) {
 return asad<4 && aget<4 && ch<64 ? int(asad*256+aget*64+ch+1) : -1;
}
void scan() {
 std::ifstream f(path,std::ios::binary);
 if(!f) throw std::runtime_error("Cannot open "+path);
 f.seekg(0,std::ios::end); const auto end=f.tellg();
 if(end<scanned) { eventFrames.clear(); scanned=12; }
 f.seekg(scanned);
 while(scanned+256<=end) {
  unsigned char h[256]; f.read(reinterpret_cast<char*>(h),256);
  Frame x{scanned,be(h+1,3)*256,be(h+8,2)*256,be(h+12,4),be(h+5,2),be(h+10,2),h[27]};
  if(h[7]!=5 || h[26]!=0 || x.asad>=4 || x.header<256 || x.size<x.header ||
     (x.type!=1 && x.type!=2) || x.item!=(x.type==1?4u:2u) ||
     uint64_t(x.items)*x.item>x.size-x.header)
   throw std::runtime_error("Invalid GET frame at byte "+std::to_string(scanned));
  if(scanned+x.size>end) break; // File may still be growing; retry on Refresh.
  eventFrames[be(h+22,4)].push_back(x);
  scanned+=x.size; f.seekg(scanned);
 }
}
void geometry(TH2Poly *hist) {
 // Same bin ordering and orientation as the legacy EventDisplay.
 for(int sectionZ=0;sectionZ<4;++sectionZ)
  for(int sectionY=0;sectionY<4;++sectionY)
   for(int col=7;col>=0;--col)
    for(int row=0;row<8;++row) {
     double z=-48+24*sectionZ+3*col, y=48-24*sectionY-3*row;
     double x[5]={z+.1,z+.1,z+2.9,z+2.9,z+.1};
     double yy[5]={y-2.9,y-.1,y-.1,y-2.9,y-2.9};
     hist->AddBin(5,x,yy);
    }
}
void setup() {
 if(canvas) return;
 charge=new TH2Poly("corea_charge","peak - baseline;z [mm];y [mm]",-60,60,-60,60);
 time=new TH2Poly("corea_time","peak time;z [mm];y [mm]",-60,60,-60,60);
 geometry(charge); geometry(time);
 pulse=new TH1D("corea_pulse","waveform;Time bucket;ADC",512,-.5,511.5);
 pulse->SetDirectory(nullptr);
 canvas=new TCanvas("corea_display","GET display",1500,520); canvas->Divide(3);
 canvas->GetPad(1)->AddExec("corea_hover_adc", "corea_hover()");
 canvas->GetPad(2)->AddExec("corea_hover_time", "corea_hover()");
}
void waveform(unsigned id) {
 if(!canvas || id<1 || id>1024) return;
 pulse->Reset();
 for(unsigned tb=0;tb<512;++tb) if(present[(id-1)*512+tb]) pulse->SetBinContent(tb+1,samples[(id-1)*512+tb]);
 pulse->SetTitle(Form("Pad %u, AsAd %u AGET %u channel %u;Time bucket;ADC",id,(id-1)/256,((id-1)%256)/64,(id-1)%64));
 canvas->cd(3); pulse->Draw("hist"); canvas->Modified(); canvas->Update();
}
bool show(uint32_t ev) {
 scan(); auto it=eventFrames.find(ev);
 if(it==eventFrames.end()) { std::cerr<<"Event "<<ev<<" not available.\n"; return false; }
 unsigned mask=0;
 for(const auto &x:it->second) { if(mask&(1u<<x.asad)) throw std::runtime_error("Duplicate AsAd frame"); mask|=1u<<x.asad; }
 if(mask!=15) { std::cerr<<"Event "<<ev<<" incomplete, AsAd mask="<<mask<<"; Refresh to retry.\n"; return false; }
 setup(); std::fill(samples.begin(),samples.end(),0); std::fill(present.begin(),present.end(),0);
 charge->Reset(""); time->Reset(""); pulse->Reset();
 std::ifstream f(path,std::ios::binary); uint64_t words=0,mapped=0;
 for(const auto &x:it->second) {
  f.seekg(x.pos+x.header); std::vector<unsigned char> raw(size_t(x.items)*x.item);
  if(!raw.empty() && !f.read(reinterpret_cast<char*>(raw.data()),raw.size())) throw std::runtime_error("Short payload");
  unsigned count[4]={0,0,0,0};
  for(unsigned i=0;i<x.items;++i) {
   uint32_t v=be(raw.data()+size_t(i)*x.item,x.item);
   unsigned a,ch,tb,adc;
   if(x.type==1) { a=(v>>30)&3; ch=(v>>23)&127; tb=(v>>14)&511; adc=v&4095; }
   else { a=(v>>14)&3; ch=count[a]%68; tb=count[a]/68; adc=v&4095; ++count[a]; }
   if(ch>=68 || tb>=512) throw std::runtime_error("Channel/time bucket out of range");
   ++words; int id=pad(x.asad,a,ch); if(id<1) continue;
   unsigned offset=(id-1)*512+tb; samples[offset]=adc; present[offset]=1; ++mapped;
  }
 }
 bestPad=1; double best=-1;
 for(unsigned id=1;id<=1024;++id) {
  double sum=0,base=0,peak=-1; unsigned nb=0,ns=0,nt=0,peakTime=0;
  for(unsigned tb=0;tb<512;++tb) {
   unsigned off=(id-1)*512+tb; if(!present[off]) continue;
   double v=samples[off]; sum+=v; ++ns;
   if(tb<100) { base+=v; ++nb; }
   if(v>peak) { peak=v; peakTime=tb; } ++nt;
  }
  if(!nt) continue;
  double height=peak-(nb?base/nb:sum/ns);
  charge->SetBinContent(id,height);
  if(height>threshold) time->SetBinContent(id,peakTime);
  if(height>best) { best=height; bestPad=id; }
 }
 current=ev;
 canvas->SetTitle(Form("run %04d raw event %u (%llu ADC words)",run,ev,(unsigned long long)words));
 canvas->cd(1); charge->Draw("colz"); canvas->cd(2); time->Draw("colz"); waveform(bestPad);
 std::cout<<"run="<<run<<" raw_event="<<ev<<" AsAd=4 ADC_words="<<words<<" mapped_samples="<<mapped<<" peak_pad="<<bestPad<<" peak_height="<<best<<std::endl;
 if(!words) std::cout<<"Empty event: no ADC samples in these frames.\n";
 return true;
}
}
// Called by ROOT when the mouse moves over either pad map.
void corea_hover() {
 using namespace CoreaOnline;
 TVirtualPad *source=gPad;
 if(!source || !canvas || source->GetCanvas()!=canvas ||
    (source!=canvas->GetPad(1) && source!=canvas->GetPad(2))) return;
 if(source->GetEvent()!=kMouseMotion && source->GetEvent()!=kButton1Down) return;
 TObject *selected=source->GetSelected();
 if(selected!=charge && selected!=CoreaOnline::time) return;
 auto *hist=static_cast<TH2Poly*>(selected);
 const double x=source->PadtoX(source->AbsPixeltoX(source->GetEventX()));
 const double y=source->PadtoY(source->AbsPixeltoY(source->GetEventY()));
 const int id=hist->FindBin(x,y);
 if(id<1 || id>1024) return;
 waveform(id);
 source->cd(); // Keep ROOT's mouse event attached to its originating pad.
}
void corea_event(unsigned event) { CoreaOnline::show(event); }
void corea_pulse(unsigned pad=0) { CoreaOnline::waveform(pad?pad:CoreaOnline::bestPad); }
void corea_next() {
 using namespace CoreaOnline; scan(); auto it=eventFrames.upper_bound(current);
 if(it!=eventFrames.end()) show(it->first); else std::cout<<"No next event yet; Refresh after new data arrive.\n";
}
void corea_previous() {
 using namespace CoreaOnline; auto it=eventFrames.lower_bound(current);
 if(it!=eventFrames.begin()) { --it; show(it->first); }
}
void corea_refresh() { CoreaOnline::show(CoreaOnline::current); }
void corea_latest() {
 using namespace CoreaOnline; scan();
 for(auto it=eventFrames.rbegin();it!=eventFrames.rend();++it) {
  unsigned mask=0; for(const auto &x:it->second) mask|=1u<<x.asad;
  if(mask==15) { show(it->first); return; }
 }
}
void online(int runNumber=24,unsigned event=0,double cut=50,const char *dataDir="/home/axis/ganacq_manip/corea/acquisition/run/cobo0") {
 using namespace CoreaOnline;
 run=runNumber; threshold=cut; path=std::string(dataDir)+Form("/run_%04d.dat",run);
 eventFrames.clear(); scanned=12; current=event; show(event);
 if(!gROOT->IsBatch() && !controls) {
  controls=new TControlBar("horizontal","Event display");
  controls->AddButton("Previous","corea_previous()"); controls->AddButton("Next","corea_next()");
  controls->AddButton("Refresh","corea_refresh()"); controls->AddButton("Latest","corea_latest()");
  controls->AddButton("Peak pulse","corea_pulse()"); controls->Show();
 }
}
