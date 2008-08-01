//---------------------------------------------------------------------------

#ifndef segment_treeH
#define segment_treeH



//int min_int(int a,int b)
//{
//  if (a<b)
//    return a;
//  else
//    return b;
//};


#ifdef BORLAND_WIN
#include <vcl.h>
#endif

class SegmentTree;

class RoundRobinInterval
{
  friend class SegmentTree;
  int modulo;
  int a, b;
protected:
  RoundRobinInterval(int m)
  {
    modulo=m;
  };
public:
  void set(int _a,int _b)
  {
    a=_a;
    b=_b;
  };


  // anf indicates the first position of the Interval,
  // end the last.
  // end<anf is allowed (Round Robin)
  // Negativ Values are handeled in a modulo Fashion

  void exclude(int anf,int end)
  {
  //XXX! Should this say a = ...; b = ...;?
  //As is, it does nothing.
       
     anf=anf % modulo;
     end=end % modulo;


  };
};



class SegmentTree
{
//  class Aktion;
//  class Interval;

  class Aktion
  {
   public:
    int typ;  //0==set, 1==add, 2=getMin
    int l;
    int r;
    int wert;
    Aktion(int _typ,int _l,int _r,int _wert)
    {
      typ=_typ;
      l=_l;
      r=_r;
      wert=_wert;
    };
    virtual ~Aktion() {};
  };

class Interval
{
  int r1; //r1>=r  [l,r1] ist das Intervall im 2^k-Baum
  public:
  int pos;//In the data array
  int l;
  int r;
  int max_r;

  int left(void)
  {
    return 2*pos+1;
  };
  int right(void)
  {
    return 2*pos+2;
  };

  virtual ~Interval() {};
  Interval(int _l, int _r, int _max_r)
  {
    pos=0;
    l=_l;
    r1=_r;
    max_r=_max_r;
    r=r1;
    if (r>max_r)
      r=max_r;
  };
  /*
  bool is_intersected_by(int a,int b)
  {
    if ( (l<=a && a<=r) || (l<=b && b<=r))
      return true;
    return false;
  };
    */
  bool is_not_intersected_by(int a,int b)
  {
    if (b<l || r<a)
      return true;
    return false;
  };

  bool interval_covered(int left,int right)
  {
    if (left<=l && right>=r)
      return(true);
    return(false);
  };

  void go_left(void)
  {
    r1=(l+r1)/2;
    r=(r1<max_r)?(r1):(max_r);
    pos=left();
  };
  void go_right(void)
  {
    l=(l+r1)/2+1;
    pos=right();
  };
  void go_up(void)
  {
    if (pos/2*2==pos) //gerade => rechtes Kind
    {
      l=l-(r1-l+1);
    }
    else
    {
      r1=r1+(r1-l+1);
      r=(r1<max_r)?(r1):(max_r);
    }

    pos=(pos-1)/2;
  };

  bool isLeaf(void)
  {
    if (l==r1)
      return true;
    return false;
  };

  bool right_child_exists(void)
  {
    if ((l+r1)/2+1<=max_r)
      return true;
    return false;
  };
};

//begin class SegmentTree

  int * data;   // Ein möglicher update-wert ist bereits enthalten
  int * update; // Dieser Wert muss nur auf die Teilbäume angewandt werden
  int * minPos;
  int size;    // Anzahl der Blätter in 2^k gequantelt
  int max_r;   // Anzahl der eingefügten Elemente -1 = Aufrufe von extend()
  int modulo;  // Anzahl der eingefügten Elemente
  int ergPos;
  int ergValue;

  void add(Interval *seg,Aktion *ak, int updateValue);

  //Interval *seg;
  //Aktion *ak;

public:

  SegmentTree(int _size)
  {
    size=1;
    while(size<_size)
      size*=2;

    data=new int[2*size];
    update=new int[2*size];
    minPos=new int[2*size];

    max_r=-1;
    modulo=0;
    //seg=new Interval(0,0,0);
    //ak=new Aktion(0,0,0,0);
  };
  ~SegmentTree()
  {
    delete []data;
    delete []update;
    delete []minPos;

    //delete seg;
    //delete ak;
  };

  #ifdef BORLAND_WIN
  AnsiString print(void);
  #endif

  int getMinPos(int l,int r); // Position des Minimums
  int getLastValue(void)
  {
    return ergValue;
  };

  void add(int l,int r,int wert);

  void set(int pos, int wert);

  void extend(int wert);

};


//---------------------------------------------------------------------------
#endif
