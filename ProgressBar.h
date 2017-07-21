// Progress bar class
// Author: Yi Chen

#include <iostream>
#include <iomanip>
#include <ostream>
#include <cstdlib>

class ProgressBar
{
private:
   std::ostream *Out;
   double Max;
   double Min;
   int Column;
   double Progress;
   int Style;
   void SanityCheck();
public:
   ProgressBar(std::ostream &out, double max = 100, double min = 0, int column = 80)
      : Out(&out), Max(max), Min(min), Column(column), Progress(0), Style(0) {SanityCheck();}
   ProgressBar(std::ostream *out, double max = 100, double min = 0, int column = 80)
      : Out(out), Max(max), Min(min), Column(column), Progress(0), Style(0) {SanityCheck();}
   ~ProgressBar() {}
   void Print();
   void PrintWithMod(int Mod);
   void Print(double progress);
   void ChangeLine() {*Out << std::endl;}
   void PrintLine() {*Out << std::endl;}
   void Update(double progress) {SetProgress(progress);}
   void Increment(double change = 1) {Progress = Progress + change;}
public:
   double GetMin() {return Min;}
   double GetMax() {return Max;}
   double GetProgress() {return Progress;}
   int GetColumn() {return Column;}
   int GetStyle() {return Style;}
   std::ostream *GetStream() {return Out;}
   double GetPercentage() {return (Progress - Min) / (Max - Min);}
public:
   void SetMin(double min) {Min = min;   SanityCheck();}
   void SetMax(double max) {Max = max;   SanityCheck();}
   void SetProgress(double progress) {Progress = progress;   SanityCheck();}
   void SetColumn(int column) {Column = column;   SanityCheck();}
   void SetStyle(int style) {if(style == -1) Style = rand() % 6; else Style = style;   SanityCheck();}
   void SetStream(std::ostream &out) {Out = &out;   SanityCheck();}
   void SetStream(std::ostream *out) {Out = out;   SanityCheck();}
};

void ProgressBar::SanityCheck()
{
   if(Min == Max)
   {
      std::cerr << "[ProgressBar] Sanity check on range failed.  Resetting to 0-100" << std::endl;
      Min = 0;
      Max = 100;
      Progress = 0;
   }
   if(Max < Min)
   {
      std::cerr << "[ProgressBar] Min > Max!  Reversing the role of the two" << std::endl;
      std::swap(Min, Max);
   }

   if(Progress < Min)
   {
      std::cerr << "[ProgressBar] Negative progress.  Resetting to minimum value" << std::endl;
      Progress = Min;
   }
   if(Progress > Max)
   {
      std::cerr << "[ProgressBar] Past-complete progress.  Resetting to maximum value" << std::endl;
      Progress = Max;
   }

   if(Column < 15)
   {
      std::cerr << "[ProgressBar] Too few columns to display the progress bar.  Set to 20" << std::endl;
      Column = 20;
   }
   if(Column > 100)
   {
      std::cerr << "[ProgressBar] Too many columns to display the progress bar.  Set to 100" << std::endl;
      Column = 100;
   }

   if(Style < 0 || Style > 7)
   {
      std::cerr << "[ProgressBar] Style invalid.  Set to a random style." << std::endl;
      std::cerr << std::endl;
      std::cerr << "FYI: available styles look like these" << std::endl;
      std::cerr << "0: [==============================>                 ]  55%" << std::endl;
      std::cerr << "1: [                            ><>                 ]  55%" << std::endl;
      std::cerr << "2: [ooooooooooooooooooooooooooooooo                 ]  55%" << std::endl;
      std::cerr << "3: [~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|                 ]  55%" << std::endl;
      std::cerr << "4: [                  <=============================]  55%" << std::endl;
      std::cerr << "5: [                  <><                           ]  55%" << std::endl;
      std::cerr << "6: Current progress: 553/1000 (55%)" << std::endl;
      std::cerr << "7: Current progress: 553" << std::endl;
      std::cerr << std::endl;

      Style = std::rand() % 8;
   }

   if(Out == NULL)
   {
      std::cerr << "[ProgressBar] Output stream is NULL.  Set to cout" << std::endl;
      Out = &std::cout;
   }
}

void ProgressBar::Print()
{
   Print(Progress);
}

void ProgressBar::PrintWithMod(int Mod)
{
   if((int)Progress % Mod == 0)
      Print(Progress);
}

void ProgressBar::Print(double progress)
{
   if(Style == 0)
   {
      int AvailableColumn = Column - 2 - 5;
      int FilledColumn = (int)(AvailableColumn * (progress - Min) / (Max - Min));

      *Out << "\033[1G[";
      for(int i = 0; i < FilledColumn - 1; i++)
         *Out << "=";
      if(FilledColumn >= 1)
         *Out << ">";
      for(int i = 0; i < AvailableColumn - FilledColumn; i++)
         *Out << " ";
      *Out << "] ";
      *Out << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5);
      *Out << "%" << std::flush;
   }
   if(Style == 1)
   {
      int AvailableColumn = Column - 2 - 5;
      int FilledColumn = (int)(AvailableColumn * (progress - Min) / (Max - Min));

      *Out << "\033[1G[";
      for(int i = 0; i < FilledColumn - 3; i++)
         *Out << " ";
      if(FilledColumn >= 3)
         *Out << ">";
      if(FilledColumn >= 2)
         *Out << "<";
      if(FilledColumn >= 1)
         *Out << ">";
      for(int i = 0; i < AvailableColumn - FilledColumn; i++)
         *Out << " ";
      *Out << "] ";
      *Out << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5);
      *Out << "%" << std::flush;
   }
   if(Style == 2)
   {
      int AvailableColumn = Column - 2 - 5;
      int FilledColumn = (int)(AvailableColumn * (progress - Min) / (Max - Min));

      *Out << "\033[1G[";
      for(int i = 0; i < FilledColumn; i++)
         *Out << "o";
      for(int i = 0; i < AvailableColumn - FilledColumn; i++)
         *Out << " ";
      *Out << "] ";
      *Out << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5);
      *Out << "%" << std::flush;
   }
   if(Style == 3)
   {
      int AvailableColumn = Column - 2 - 5;
      int FilledColumn = (int)(AvailableColumn * (progress - Min) / (Max - Min));

      *Out << "\033[1G[";
      for(int i = 0; i < FilledColumn - 1; i++)
         *Out << "~";
      if(FilledColumn >= 1)
         *Out << "|";
      for(int i = 0; i < AvailableColumn - FilledColumn; i++)
         *Out << " ";
      *Out << "] ";
      *Out << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5);
      *Out << "%" << std::flush;
   }
   if(Style == 4)
   {
      int AvailableColumn = Column - 2 - 5;
      int FilledColumn = (int)(AvailableColumn * (progress - Min) / (Max - Min));

      *Out << "\033[1G[";
      for(int i = 0; i < AvailableColumn - FilledColumn; i++)
         *Out << " ";
      if(FilledColumn >= 1)
         *Out << "<";
      for(int i = 0; i < FilledColumn - 1; i++)
         *Out << "=";
      *Out << "] ";
      *Out << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5);
      *Out << "%" << std::flush;
   }
   if(Style == 5)
   {
      int AvailableColumn = Column - 2 - 5;
      int FilledColumn = (int)(AvailableColumn * (progress - Min) / (Max - Min));

      *Out << "\033[1G[";
      for(int i = 0; i < AvailableColumn - FilledColumn; i++)
         *Out << " ";
      if(FilledColumn >= 1)
         *Out << "<";
      if(FilledColumn >= 2)
         *Out << ">";
      if(FilledColumn >= 3)
         *Out << "<";
      for(int i = 0; i < FilledColumn - 3; i++)
         *Out << " ";
      *Out << "] ";
      *Out << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5);
      *Out << "%" << std::flush;
   }
   if(Style == 6)
   {
      *Out << "\033[1GCurrent progress: " << progress - Min << "/" << Max - Min << " ("
         << std::setw(3) << std::setfill(' ') << (int)((progress - Min) / (Max - Min) * 100 + 0.5)
         << "%)" << std::flush;
   }
   if(Style == 7)
      *Out << "\033[1GCurrent progress: " << progress - Min << std::flush;
}






