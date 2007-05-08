
#include <winstd.H>

#include <map>
#include <algorithm>
#include <iostream>

#include <BoxArray.H>
#include <BoxList.H>
#include <Profiler.H>

void
BoxList::join (const BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    std::list<Box> lb = blist.lbox;
    lbox.splice(lbox.end(), lb);
}

void
BoxList::catenate (BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    lbox.splice(lbox.end(), blist.lbox);
    BL_ASSERT(blist.isEmpty());
}

BoxList&
BoxList::remove (const Box& bx)
{
    BL_ASSERT(ixType() == bx.ixType());
    lbox.remove(bx);
    return *this;
}

BoxList&
BoxList::remove (iterator bli)
{
    BL_ASSERT(ixType() == bli->ixType());
    lbox.erase(bli);
    return *this;
}

BoxList
BoxLib::intersect (const BoxList& bl,
		   const Box&     b)
{
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newbl(bl);
    return newbl.intersect(b);
}

BoxList
BoxLib::intersect (const BoxList& bl,
                   const BoxList& br)
{
    BL_ASSERT(bl.ixType() == br.ixType());
    BoxList newbl(bl);
    return newbl.intersect(br);
}

BoxList
BoxLib::refine (const BoxList& bl,
		int            ratio)
{
    BoxList nbl(bl);
    return nbl.refine(ratio);
}

BoxList
BoxLib::coarsen (const BoxList& bl,
                 int            ratio)
{
    BoxList nbl(bl);
    return nbl.coarsen(ratio);
}

BoxList
BoxLib::accrete (const BoxList& bl,
                 int            sz)
{
    BoxList nbl(bl);
    return nbl.accrete(sz);
}

BoxList
BoxLib::removeOverlap (const BoxList& bl)
{
    BoxArray ba(bl);
    return ba.removeOverlap();
}

bool
BoxList::operator!= (const BoxList& rhs) const
{
    return !operator==(rhs);
}

BoxList::BoxList ()
    :
    lbox(),
    btype(IndexType::TheCellType())
{}

BoxList::BoxList (const Box& bx)
    : btype(bx.ixType())
{
    push_back(bx);
}

BoxList::BoxList (IndexType _btype)
    :
    lbox(),
    btype(_btype)
{}

BoxList::BoxList (const BoxArray &ba)
    :
    lbox(),
    btype()
{
    if (ba.size() > 0)
        btype = ba[0].ixType();
    for (int i = 0; i < ba.size(); ++i)
        push_back(ba[i]);
}

bool
BoxList::ok () const
{
    const_iterator bli = begin();
    if ( bli != end() )
    {
        for (Box b(*bli); bli != end(); ++bli)
            if (!(bli->ok() && bli->sameType(b)))
                return false;
    }
    return true;
}

bool
BoxList::isDisjoint () const
{
    for (const_iterator bli = begin(); bli != end(); ++bli)
    {
        const_iterator bli2 = bli;
        //
        // Skip the first element.
        //
        ++bli2;
        for (; bli2 != end(); ++bli2)
            if (bli->intersects(*bli2))
                return false;
    }
    return true;
}

bool
BoxList::contains (const IntVect& v) const
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::contains(IntVect)");

    for (const_iterator bli = begin(); bli != end(); ++bli)
        if (bli->contains(v))
            return true;

    return false;
}

bool
BoxList::contains (const Box& b) const
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::contains(Box)");

    if (isEmpty()) return false;

    BL_ASSERT(ixType() == b.ixType());

    BoxList bnew = BoxLib::complementIn(b,*this);

    return bnew.isEmpty();
}

bool
BoxList::contains (const BoxList&  bl) const
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::contains(BoxList)");

    if (isEmpty() || bl.isEmpty()) return false;

    BL_ASSERT(ixType() == bl.ixType());

    if (!minimalBox().contains(bl.minimalBox())) return false;

    BoxArray ba(*this);

    for (const_iterator bli = bl.begin(); bli != bl.end(); ++bli)
        if (!ba.contains(*bli))
            return false;

    return true;
}

bool
BoxList::contains (const BoxArray&  ba) const
{
    BoxArray tba(*this);
    return ba.contains(tba);
}

BoxList&
BoxList::intersect (const Box& b)
{
    BL_ASSERT(ixType() == b.ixType());
    for (iterator bli = begin(); bli != end(); )
    {
        Box bx = *bli & b;

        if (bx.ok())
        {
            *bli = bx;
            ++bli;
        }
        else
        {
            bli = lbox.erase(bli);
        }
    }
    return *this;
}

BoxList&
BoxList::intersect (const BoxList& b)
{
    BL_ASSERT(ixType() == b.ixType());

    BoxList bl(b.ixType());

    for (iterator lhs = begin(); lhs != end(); ++lhs)
    {
        for (const_iterator rhs = b.begin(); rhs != b.end(); ++rhs)
        {
            Box bx = *lhs & *rhs;
            if (bx.ok())
                bl.push_back(bx);
        }
    }

    *this = bl;

    return *this;
}

BoxList
BoxLib::complementIn (const Box&     b,
                      const BoxList& bl)
{
    BL_PROFILE("BoxLib::complementIn(Box,BoxList)");
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newb(b.ixType());
    newb.complementIn(b,bl);
    return newb;
}

BoxList&
BoxList::complementIn (const Box&     b,
                       const BoxList& bl)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::complementIn()");

    BL_ASSERT(bl.ixType() == b.ixType());

    if (bl.size() == 1)
    {
        *this = BoxLib::boxDiff(b,bl.front());
    }
    else
    {
        clear();

        Box     minbox = bl.minimalBox();
        BoxList tmpbl  = BoxLib::boxDiff(b,minbox);

        catenate(tmpbl);

        BoxArray ba(bl);

        BoxList mesh(b.ixType());
        if (minbox.ok())
            mesh.push_back(minbox);
        mesh.maxSize(BL_SPACEDIM == 3 ? 64 : 128);

        for (BoxList::const_iterator bli = mesh.begin(); bli != mesh.end(); ++bli)
        {
            const Box bx = *bli & b;

            if (!bx.ok()) continue;

            std::vector< std::pair<int,Box> > isects = ba.intersections(bx);

            if (!isects.empty())
            {
                tmpbl.clear();
                BoxList tm(b.ixType());
                for (int i = 0; i < isects.size(); i++)
                    tmpbl.push_back(isects[i].second);
                tm.complementIn_base(bx,tmpbl);
                catenate(tm);
            }
            else
            {
                push_back(bx);
            }
        }
    }

    return *this;
}

BoxList&
BoxList::complementIn_base (const Box&     b,
                            const BoxList& bl)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::complementIn_base()");

    BL_ASSERT(bl.ixType() == b.ixType());

    clear();

    push_back(b);

    for (const_iterator bli = bl.begin(); bli != bl.end() && isNotEmpty(); ++bli)
    {
        for (iterator newbli = lbox.begin(); newbli != lbox.end(); )
        {
            if (newbli->intersects(*bli))
            {
                BoxList tm = BoxLib::boxDiff(*newbli, *bli);
                lbox.splice(lbox.begin(), tm.lbox);
                lbox.erase(newbli++);
            }
            else
            {
                ++newbli;
            }
        }
    }

    return *this;
}

BoxList&
BoxList::refine (int ratio)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::refine (const IntVect& ratio)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (int ratio)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (const IntVect& ratio)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::accrete (int sz)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->grow(sz);
    }
    return *this;
}

BoxList&
BoxList::accrete (IntVect sz)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->grow(sz);
    }
    return *this;
}

BoxList&
BoxList::shift (int dir,
                int nzones)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->shift(dir, nzones);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (int dir,
                    int num_halfs)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->shiftHalf(dir, num_halfs);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (const IntVect& iv)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->shiftHalf(iv);
    }
    return *this;
}

//
// Returns a list of boxes defining the compliment of b2 in b1in.
//

BoxList
BoxLib::boxDiff (const Box& b1in,
		 const Box& b2)
{
   BL_ASSERT(b1in.sameType(b2));
  
   Box b1(b1in);
   BoxList b_list(b1.ixType());

   if ( !b2.contains(b1) )
   {
       if ( !b1.intersects(b2) )
       {
           b_list.push_back(b1);
       }
       else
       {
           const int* b2lo = b2.loVect();
           const int* b2hi = b2.hiVect();

           for (int i = 0; i < BL_SPACEDIM; i++)
           {
               const int* b1lo = b1.loVect();
               const int* b1hi = b1.hiVect();

               if ((b1lo[i] < b2lo[i]) && (b2lo[i] <= b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b1lo[i]);
                   bn.setBig(i,b2lo[i]-1);
                   b_list.push_back(bn);
                   b1.setSmall(i,b2lo[i]);
               }
               if ((b1lo[i] <= b2hi[i]) && (b2hi[i] < b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b2hi[i]+1);
                   bn.setBig(i,b1hi[i]);
                   b_list.push_back(bn);
                   b1.setBig(i,b2hi[i]);
               }
           }
       }
   }
   return b_list;
}

int
BoxList::simplify ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::simplify()");

    int count = 0;

    std::multimap<IntVect, Box, IntVect::Compare> mmap;

    typedef std::multimap<IntVect, Box, IntVect::Compare>::iterator MMapIter;

    for (iterator it = begin(); it != end(); ++it)
    {
        mmap.insert(std::make_pair<IntVect,Box>(it->smallEnd(),*it));
    }

    clear();

    for (MMapIter it = mmap.begin(); it != mmap.end(); ++it)
    {
        push_back(it->second);
    }

    const int N = 50;

    BoxList tbl(ixType());

    while (!lbox.empty())
    {
        BoxList tmp(ixType());

        while (tmp.size() < N && !lbox.empty())
        {
            tmp.push_back(lbox.front());
            lbox.pop_front();
        }

        count += tmp.simplify_doit();

        tbl.catenate(tmp);
    }

    lbox.swap(tbl.lbox);

    return count;
}

int
BoxList::simplify_doit ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::simplify_doit()");
    //
    // Try to merge adjacent boxes.
    //
    int count = 0;
    int lo[BL_SPACEDIM];
    int hi[BL_SPACEDIM];

    for (iterator bla = begin(); bla != end(); )
    {
        const int* alo   = bla->loVect();
        const int* ahi   = bla->hiVect();
        bool       match = false;
        iterator blb = bla;
        ++blb;
        while ( blb != end() )
        {
            const int* blo = blb->loVect();
            const int* bhi = blb->hiVect();
            //
            // Determine of a and b can be coalasced.
            // They must have equal extents in all index direciton
            // except possibly one and must abutt in that direction.
            //
            bool canjoin = true;
            int  joincnt = 0;
            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                if (alo[i]==blo[i] && ahi[i]==bhi[i])
                {
                    lo[i] = alo[i];
                    hi[i] = ahi[i];
                }
                else if (alo[i]<=blo[i] && blo[i]<=ahi[i]+1)
                {
                    lo[i] = alo[i];
                    hi[i] = std::max(ahi[i],bhi[i]);
                    joincnt++;
                }
                else if (blo[i]<=alo[i] && alo[i]<=bhi[i]+1)
                {
                    lo[i] = blo[i];
                    hi[i] = std::max(ahi[i],bhi[i]);
                    joincnt++;
                }
                else
                {
                    canjoin = false;
                    break;
                }
            }
            if (canjoin && (joincnt <= 1))
            {
                //
                // Modify b and remove a from the list.
                //
                blb->setSmall(IntVect(lo));
                blb->setBig(IntVect(hi));
                lbox.erase(bla++);
                count++;
                match = true;
                break;
            }
            else
            {
                //
                // No match found, try next element.
                //
                ++blb;
            }
        }
        //
        // If a match was found, a was already advanced in the list.
        //
        if (!match)
            ++bla;
    }
    return count;
}

int
BoxList::minimize ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::minimize()");
    int cnt = 0;
    for (int n; (n=simplify()) > 0; )
        cnt += n;
    return cnt;
}

Box
BoxList::minimalBox () const
{
    Box minbox(IntVect::TheUnitVector(), IntVect::TheZeroVector(), ixType());
    if ( !isEmpty() )
    {
        const_iterator bli = begin();
        minbox = *bli;
        while ( bli != end() )
	{
            minbox.minBox(*bli++);
	}
    }
    return minbox;
}

BoxList&
BoxList::maxSize (const IntVect& chunk)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::maxSize()");

    for (iterator bli = begin(); bli != end(); ++bli)
    {
        IntVect boxlen = bli->size();
        const int* len = boxlen.getVect();

        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            if (len[i] > chunk[i])
            {
                //
                // Reduce by powers of 2.
                //
                int ratio = 1;
                int bs    = chunk[i];
                int nlen  = len[i];
                while ((bs%2 == 0) && (nlen%2 == 0))
                {
                    ratio *= 2;
                    bs    /= 2;
                    nlen  /= 2;
                }
                //
                // Determine number and size of (coarsened) cuts.
                //
                const int numblk = nlen/bs + (nlen%bs ? 1 : 0);
                const int size   = nlen/numblk;
                const int extra  = nlen%numblk;
                //
                // Number of cuts = number of blocks - 1.
                //
                for (int k = 0; k < numblk-1; k++)
                {
                    //
                    // Compute size of this chunk, expand by power of 2.
                    //
                    const int ksize = (k < extra ? size+1 : size) * ratio;
                    //
                    // Chop from high end.
                    //
                    const int pos = bli->bigEnd(i) - ksize + 1;

                    push_back(bli->chop(i,pos));
                }
            }
        }
        //
        // b has been chopped down to size and pieces split off
        // have been added to the end of the list so that they
        // can be checked for splitting (in other directions) later.
        //
    }
    return *this;
}

BoxList&
BoxList::maxSize (int chunk)
{
    return maxSize(IntVect(D_DECL(chunk,chunk,chunk)));
}

BoxList&
BoxList::surroundingNodes ()
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->surroundingNodes();
    }
    return *this;
}

BoxList&
BoxList::surroundingNodes (int dir)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->surroundingNodes(dir);
    }
    return *this;
}

BoxList&
BoxList::enclosedCells ()
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->enclosedCells();
    }
    return *this;
}

BoxList&
BoxList::enclosedCells (int dir)
{
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->enclosedCells(dir);
    }
    return *this;
}

BoxList&
BoxList::convert (IndexType typ)
{
    btype = typ;
    for (iterator bli = begin(); bli != end(); ++bli)
    {
        bli->convert(typ);
    }
    return *this;
}

std::ostream&
operator<< (std::ostream&  os,
            const BoxList& blist)
{
    BoxList::const_iterator bli = blist.begin();
    os << "(BoxList " << blist.size() << ' ' << blist.ixType() << '\n';
    for (int count = 1; bli != blist.end(); ++bli, ++count)
    {
        os << count << " : " << *bli << '\n';
    }
    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,BoxList&) failed");

    return os;
}

bool
BoxList::operator== (const BoxList& rhs) const
{
    if ( !(size() == rhs.size()) ) return false;

    BoxList::const_iterator liter = begin(), riter = rhs.begin();
    for (; liter != end(); ++liter, ++riter)
        if ( !( *liter == *riter) )
            return false;
    return true;
}
