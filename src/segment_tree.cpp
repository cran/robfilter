//---------------------------------------------------------------------------

//#pragma hdrstop

#include "segment_tree.h"

//#ifndef BORLAND_WIN
//#endif

#define SEGTREE_SET 0
#define SEGTREE_ADD 1
#define SEGTREE_GETMIN 2





  #ifdef BORLAND_WIN
  AnsiString SegmentTree::print(void)
  {
    AnsiString s;
    s="(data,update,minPos):\r\n";

    int len=1;
    int pos=0;

    while(pos<2*size-1)
    {
      for(int i=0;i<len && pos<2*size-1;i++)
      {
        s=s+"("+data[pos]+","+update[pos]+","+minPos[pos]+")  ";
        pos++;
      }
      s+="\r\n";
      len*=2;
    }
    return s;
  };
  #endif

  int SegmentTree::getMinPos(int l,int r) // Position des Minimums
  {
    while(l<0)
      l+=modulo;
    while(r<0)
      r+=modulo;

    l=l%modulo;
    r=r%modulo;

    ergPos=-1;
    ergValue=0;

    if (l<=r)
    {
      Interval seg(0,size-1,max_r);
      Aktion ak(2,l,r,0);
      add(&seg,&ak,0);
    }
    else
    {
      int p1,e1,e2;

      {
       Interval seg(0,size-1,max_r);
       Aktion ak(2,l,max_r,0);
       add(&seg,&ak,0);
      }
      p1=ergPos;
      e1=ergValue;
      ergPos=-1;
      ergValue=0;
      {
       Interval seg(0,size-1,max_r);
       Aktion ak(2,0,r,0);
       add(&seg,&ak,0);
      }
      //p2=ergPos;
      e2=ergValue;

      if (e1<e2)
      {
        ergPos=p1;
        ergValue=e1;
      }
    }
    return ergPos-size+1;
  };

  void SegmentTree::add(int l,int r,int wert)
  {
    while(l<0)
      l+=modulo;
    while(r<0)
      r+=modulo;

    l=l%modulo;
    r=r%modulo;

    if (l<=r)
    {
      Interval seg(0,size-1,max_r);
      Aktion ak(1,l,r,wert);
      add(&seg,&ak,0);
    }
    else
    {
      {
       Interval seg(0,size-1,max_r);
       Aktion ak(1,l,max_r,wert);
       add(&seg,&ak,0);
      }
      {
       Interval seg(0,size-1,max_r);
       Aktion ak(1,0,r,wert);
       add(&seg,&ak,0);
      }
    }
  };

  void SegmentTree::add(Interval *seg,Aktion *ak, int updateValue)
  {
    update[seg->pos]+=updateValue;
    data[seg->pos]+=updateValue;

    if (seg->is_not_intersected_by(ak->l,ak->r))
      return;

    bool covered=seg->interval_covered(ak->l,ak->r);
    bool branch=true;

    if (covered)
    {
      branch=false;
      if (ak->typ==2)
        if (ergPos==-1 || data[seg->pos]<ergValue)
        {
          ergPos=minPos[seg->pos];
          ergValue=data[seg->pos];
        }
      if (ak->typ==1)
      {
        update[seg->pos]+=ak->wert;
        data[seg->pos]+=ak->wert;
      }
      if (ak->typ==0)
      {
        update[seg->pos]=0;
        minPos[seg->pos]=seg->pos;

        if (seg->isLeaf())
          data[seg->pos]=ak->wert;
        else
          branch=true;
      }
    }

    if (branch)
    {
      updateValue=update[seg->pos];
      update[seg->pos]=0;

      seg->go_left();
      add(seg,ak,updateValue);
      seg->go_up();

      if (seg->right_child_exists())
      {
        seg->go_right();
        add(seg,ak,updateValue);
        seg->go_up();
      }
      // Werte von unten einsammeln
      if (!seg->isLeaf())
      {
        if (!seg->right_child_exists() || data[seg->left()]<data[seg->right()])
        {
          data[seg->pos]=data[seg->left()];
          minPos[seg->pos]=minPos[seg->left()];
        }
        else
        {
          data[seg->pos]=data[seg->right()];
          minPos[seg->pos]=minPos[seg->right()];
        }
      }
    }

    return;
  };

  void SegmentTree::set(int pos, int wert)
  {
    Interval seg(0,size-1,max_r);
    Aktion ak(0,pos,pos,wert);

    add(&seg,&ak,0);
  };

  void SegmentTree::extend(int wert)
  {
    max_r++;
    modulo++;
    set(max_r,wert);
  };


//---------------------------------------------------------------------------

//#pragma package(smart_init)
