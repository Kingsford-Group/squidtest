#include "SV.h"
#include "GtfTrans.h"

using namespace std;

SV_t::SV_t(map<string,int> & RefTable, int numfile, string files[]){
    updated=false;
    vSimpleSV.reserve(2048);
    int count=-1;
    for(int i=0;i<numfile;i++){
        string line;
        ifstream input(files[i]);
        getline(input,line);
        vector<string> strs;
        vector<string> namestrs;
        boost::split(strs,line,boost::is_any_of("\t"));
        if(strs.size()==6){ //Deletions
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                boost::split(namestrs, strs[1], boost::is_any_of(" "));
                SimpleSV_t tmp(RefTable[namestrs[0]], stoi(strs[2])-1, stoi(strs[3]), DEL, ++count);
                vSimpleSV.push_back(tmp);
                if(vSimpleSV.capacity()==vSimpleSV.size())
                    vSimpleSV.reserve(2*vSimpleSV.size());
            }
        }
        else if(strs.size()==7 && strs[5]=="Duplications"){ //Duplications, the first copy is normal not duplication, start from the second copy
            while(getline(input,line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                boost::split(namestrs, strs[1], boost::is_any_of(" "));
                if(strs[5]=="1")
                    continue;
                SimpleSV_t tmp(RefTable[namestrs[0]], stoi(strs[3]), (stoi(strs[3])-stoi(strs[2])+1)*(stoi(strs[5])-1), INS, ++count);
                vSimpleSV.push_back(tmp);
                if(vSimpleSV.capacity()==vSimpleSV.size())
                    vSimpleSV.reserve(2*vSimpleSV.size());
            }
        }
        else if(strs.size()==7){ //Inversions
            while(getline(input,line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                boost::split(namestrs, strs[1], boost::is_any_of(" "));
                SimpleSV_t tmp(RefTable[namestrs[0]], stoi(strs[2])-1, stoi(strs[3]), INV, ++count);
                vSimpleSV.push_back(tmp);
                if(vSimpleSV.capacity()==vSimpleSV.size())
                    vSimpleSV.reserve(2*vSimpleSV.size());
            }
        }
        else if(strs.size()==12 && strs[8]=="Copied"){ //Insertions, which cuts or copies from chrA and paste to chr B, aka transposon. =deletion+insertion
            while(getline(input,line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                boost::split(namestrs, strs[4], boost::is_any_of(" "));
                SimpleSV_t tmp1(RefTable[namestrs[0]], stoi(strs[5])-1, stoi(strs[6])-stoi(strs[5])+1, INS, ++count);
                vSimpleSV.push_back(tmp1);
                if(vSimpleSV.capacity()==vSimpleSV.size())
                    vSimpleSV.reserve(2*vSimpleSV.size());
                if(boost::to_upper_copy<std::string>(strs[8])=="FALSE"){
                    boost::split(namestrs, strs[1], boost::is_any_of(" "));
                    SimpleSV_t tmp2(RefTable[namestrs[0]], stoi(strs[2])-1, stoi(strs[3]), DEL, count);
                    vSimpleSV.push_back(tmp2);
                    if(vSimpleSV.capacity()==vSimpleSV.size())
                        vSimpleSV.reserve(2*vSimpleSV.size());
                }
            }
        }
        else if(strs.size()==12 && strs[9]=="Balanced"){
            while(getline(input,line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                string chr1, chr2;
                boost::split(namestrs, strs[1], boost::is_any_of(" ")); chr1=namestrs[0];
                boost::split(namestrs, strs[5], boost::is_any_of(" ")); chr2=namestrs[0];
                int mid1=stoi(strs[2]), mid2=stoi(strs[6]);
                int Pos1=(mid1==1)?stoi(strs[3]):(mid1-1), Pos2=(mid2==1)?stoi(strs[7]):(mid2-1);
                DirType DT1=(mid1==1)?(DirType::right):(DirType::left), DT2=(mid2==1)?(DirType::right):(DirType::left);
                TRA_t tmp(RefTable[chr1], Pos1, DT1, RefTable[chr2], Pos2, DT2, ++count);
                tmp.DecideSwap();
                vTRA.push_back(tmp);
                if(vTRA.capacity()==vTRA.size())
                    vTRA.reserve(2*vTRA.size());
            }
        }
        input.close();
    }
    vSimpleSV.reserve(vSimpleSV.size()); vTRA.reserve(vTRA.size());
    num=count+1;
};

void SV_t::Update2NewPos(map<int,int>& RefLength){
    //Using the fact that RSVSim only simulate non-overlap SVs
    vector<TRA_t>::iterator ittra;
    vector<SimpleSV_t>::iterator itsim;
    for(ittra=vTRA.begin(); ittra!=vTRA.end(); ittra++){
        for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
            ittra->UpdateSimpleSV(RefLength, itsim);
        }
        for(vector<TRA_t>::iterator tmpit=vTRA.begin(); tmpit!=vTRA.end(); tmpit++){
            if(ittra==tmpit)
                continue;
            ittra->UpdateTRA(RefLength, tmpit);
        }
        ittra->EditnReverse(RefLength);
    }
    for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
        for(vector<SimpleSV_t>::iterator tmpit=vSimpleSV.begin(); tmpit!=vSimpleSV.end(); tmpit++){
            if(itsim==tmpit)
                continue;
            itsim->UpdateSimpleSV(tmpit);
        }
        for(ittra=vTRA.begin(); ittra!=vTRA.end(); ittra++){
            itsim->UpdateTRA(ittra);
        }
        itsim->EditnReverse(RefLength);
    }
    updated=!updated;
};

vector< tuple<int,int,int,int> > SV_t::UpdateNextRound2new(SV_t svnext, map<int,int>& RefLength){
    vector< tuple<int,int,int,int> > BreakPoint; BreakPoint.reserve(2*(vSimpleSV.size()+vTRA.size()+svnext.vSimpleSV.size()+svnext.vTRA.size()));
    int IDoffset=num;
    if(!updated)
        Update2NewPos(RefLength);
    if(svnext.updated)
        svnext.Update2NewPos(RefLength);
    vector<TRA_t>::iterator ittra;
    vector<SimpleSV_t>::iterator itsim;
    // Update current SV by next round SV cannot assume non-overlapping, in fact they are overlapping
    // For overlapped SV, only care about the breakpoint.
    for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
        BreakPoint.push_back(make_tuple(itsim->RefID, itsim->StartPos, itsim->ID, (int)itsim->Type));
        if(itsim->Type!=INS)
            BreakPoint.push_back(make_tuple(itsim->RefID, itsim->EndPos, itsim->ID, (int)itsim->Type));
    }
    for(ittra=vTRA.begin(); ittra!=vTRA.end(); ittra++){
        BreakPoint.push_back(make_tuple(ittra->Ref1, ittra->Pos1, ittra->ID, 3));
        BreakPoint.push_back(make_tuple(ittra->Ref2, ittra->Pos2, ittra->ID, 3));
    }
    for(itsim=svnext.vSimpleSV.begin(); itsim!=svnext.vSimpleSV.end(); itsim++){
        for(vector< tuple<int,int,int,int> >::iterator it=BreakPoint.begin(); it!=BreakPoint.end(); it++){
            pair<int,int> newBP=itsim->UpdatePoint(make_pair(get<0>(*it),get<1>(*it)));
            get<0>(*it)=newBP.first; get<1>(*it)=newBP.second;
        }
        for(vector<SimpleSV_t>::iterator tmpit=svnext.vSimpleSV.begin(); tmpit!=svnext.vSimpleSV.end(); tmpit++){
            if(itsim==tmpit)
                continue;
            itsim->UpdateSimpleSV(tmpit);
        }
        for(ittra=svnext.vTRA.begin(); ittra!=svnext.vTRA.end(); ittra++){
            itsim->UpdateTRA(ittra);
        }
        itsim->EditnReverse(RefLength);
    }
    for(ittra=svnext.vTRA.begin(); ittra!=svnext.vTRA.end(); ittra++){
        for(vector< tuple<int,int,int,int> >::iterator it=BreakPoint.begin(); it!=BreakPoint.end(); it++){
            pair<int,int> newBP=ittra->UpdatePoint(RefLength, make_pair(get<0>(*it),get<1>(*it)), DirType::left);
            get<0>(*it)=newBP.first; get<1>(*it)=newBP.second;
        }
        for(itsim=svnext.vSimpleSV.begin(); itsim!=svnext.vSimpleSV.end(); itsim++){
            ittra->UpdateSimpleSV(RefLength, itsim);
        }
        for(vector<TRA_t>::iterator tmpit=svnext.vTRA.begin(); tmpit!=svnext.vTRA.end(); tmpit++){
            if(ittra==tmpit)
                continue;
            ittra->UpdateTRA(RefLength, tmpit);
        }
        ittra->EditnReverse(RefLength);
    }
    for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
        if(itsim->Type==INS)
            BreakPoint.push_back(make_tuple(itsim->RefID, itsim->EndPos, itsim->ID, (int)itsim->Type));
    }
    for(itsim=svnext.vSimpleSV.begin(); itsim!=svnext.vSimpleSV.end(); itsim++){
        BreakPoint.push_back(make_tuple(itsim->RefID, itsim->StartPos, IDoffset+itsim->ID, (int)itsim->Type));
        BreakPoint.push_back(make_tuple(itsim->RefID, itsim->EndPos, IDoffset+itsim->ID, (int)itsim->Type));
    }
    for(ittra=svnext.vTRA.begin(); ittra!=svnext.vTRA.end(); ittra++){
        BreakPoint.push_back(make_tuple(ittra->Ref1, ittra->Pos1, IDoffset+ittra->ID, 3));
        BreakPoint.push_back(make_tuple(ittra->Ref2, ittra->Pos2, IDoffset+ittra->ID, 3));
    }
    sort(BreakPoint.begin(), BreakPoint.end(), SV_t::tupleCompare);
    vector< tuple<int,int,int,int> >::iterator it=unique(BreakPoint.begin(), BreakPoint.end(), SV_t::tupleEqual);
    BreakPoint.resize(distance(BreakPoint.begin(), it));
    return BreakPoint;
};

vector< tuple<int,int,int,int> > SV_t::UpdateNextRound2old(SV_t svnext, map<int,int>& RefLength){
    vector< tuple<int,int,int,int> > BreakPoint; BreakPoint.reserve(2*(vSimpleSV.size()+vTRA.size()+svnext.vSimpleSV.size()+svnext.vTRA.size()));
    int IDoffset=num;
    if(!updated)
        Update2NewPos(RefLength);
    if(svnext.updated)
        svnext.Update2NewPos(RefLength);
    vector<SimpleSV_t>::iterator itsim;
    vector<TRA_t>::iterator ittra;
    for(itsim=svnext.vSimpleSV.begin(); itsim!=svnext.vSimpleSV.end(); itsim++){
        BreakPoint.push_back(make_tuple(itsim->RefID, itsim->StartPos, IDoffset+itsim->ID, (int)itsim->Type));
        if(itsim->Type!=INS)
            BreakPoint.push_back(make_tuple(itsim->RefID, itsim->EndPos, IDoffset+itsim->ID, (int)itsim->Type));
    }
    for(ittra=svnext.vTRA.begin(); ittra!=svnext.vTRA.end(); ittra++){
        BreakPoint.push_back(make_tuple(ittra->Ref1, ittra->Pos1, IDoffset+ittra->ID, 3));
        BreakPoint.push_back(make_tuple(ittra->Ref2, ittra->Pos2, IDoffset+ittra->ID, 3));
    }
    for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
        for(vector< tuple<int,int,int,int> >::iterator it=BreakPoint.begin(); it!=BreakPoint.end(); it++){
            pair<int,int> newBP=itsim->UpdatePoint(make_pair(get<0>(*it),get<1>(*it)));
            get<0>(*it)=newBP.first; get<1>(*it)=newBP.second;
        }
        for(vector<SimpleSV_t>::iterator tmpit=vSimpleSV.begin(); tmpit!=vSimpleSV.end(); tmpit++){
            if(itsim==tmpit)
                continue;
            itsim->UpdateSimpleSV(tmpit);
        }
        for(ittra=vTRA.begin(); ittra!=vTRA.end(); ittra++){
            itsim->UpdateTRA(ittra);
        }
        itsim->EditnReverse(RefLength);
    }
    for(ittra=vTRA.begin(); ittra!=vTRA.end(); ittra++){
        for(vector< tuple<int,int,int,int> >::iterator it=BreakPoint.begin(); it!=BreakPoint.end(); it++){
            pair<int,int> newBP=ittra->UpdatePoint(RefLength, make_pair(get<0>(*it),get<1>(*it)), DirType::left);
            get<0>(*it)=newBP.first; get<1>(*it)=newBP.second;
        }
        for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
            ittra->UpdateSimpleSV(RefLength, itsim);
        }
        for(vector<TRA_t>::iterator tmpit=vTRA.begin(); tmpit!=vTRA.end(); tmpit++){
            if(ittra==tmpit)
                continue;
            ittra->UpdateTRA(RefLength, tmpit);
        }
        ittra->EditnReverse(RefLength);
    }
    for(itsim=svnext.vSimpleSV.begin(); itsim!=svnext.vSimpleSV.end(); itsim++){
        if(itsim->Type==INS)
            BreakPoint.push_back(make_tuple(itsim->RefID, itsim->EndPos, IDoffset+itsim->ID, (int)itsim->Type));
    }
    for(itsim=vSimpleSV.begin(); itsim!=vSimpleSV.end(); itsim++){
        BreakPoint.push_back(make_tuple(itsim->RefID, itsim->StartPos, itsim->ID, (int)itsim->Type));
        BreakPoint.push_back(make_tuple(itsim->RefID, itsim->EndPos, itsim->ID, (int)itsim->Type));
    }
    for(ittra=vTRA.begin(); ittra!=vTRA.end(); ittra++){
        BreakPoint.push_back(make_tuple(ittra->Ref1, ittra->Pos1, ittra->ID, 3));
        BreakPoint.push_back(make_tuple(ittra->Ref2, ittra->Pos2, ittra->ID, 3));
    }
    sort(BreakPoint.begin(), BreakPoint.end(), tupleCompare);
    vector< tuple<int,int,int,int> >::iterator endit=unique(BreakPoint.begin(), BreakPoint.end(), tupleEqual);
    BreakPoint.resize(distance(BreakPoint.begin(), endit));
    return BreakPoint;
};

void SV_t::WriteNextRoundBPPos(vector< tuple<int,int,int,int> >& BreakPoint, char* outputfile){
    int sizeBP=(int) BreakPoint.size();
    FILE * fp=fopen(outputfile, "wb");
    if(fp!=NULL){
        fwrite(&sizeBP, sizeof(int), 1, fp);
        for(int i=0; i<sizeBP; i++){
            fwrite(&(get<0>(BreakPoint[i])), sizeof(int), 1, fp);
            fwrite(&(get<1>(BreakPoint[i])), sizeof(int), 1, fp);
            fwrite(&(get<2>(BreakPoint[i])), sizeof(int), 1, fp);
            //fwrite(&(get<3>(BreakPoint[i])), sizeof(int), 1, fp);
        }
        fclose(fp);
    }
    else
        cout<<"Cannot open output file \""<<outputfile<<"\"\n";
};

void SV_t::WriteNextRoundBPPos(vector<string>& RefName, vector< tuple<int,int,int,int> >& BreakPoint, char* outputfile){
    int sizeBP=(int) BreakPoint.size();
    ofstream output(outputfile, ios::out);
    if(output.is_open()){
        cout<<BreakPoint.size()<<endl;
        output<<"##chr\tpos1\tpos2\tID\ttype\n";
        int curid=-1, i=0, curtype=-1;
        do{
            curid=get<2>(BreakPoint[i]); curtype=get<3>(BreakPoint[i]);
            output<<RefName[get<0>(BreakPoint[i])]<<'\t'<<get<1>(BreakPoint[i]);
            i++;
            if(i<sizeBP && get<2>(BreakPoint[i])==curid && get<3>(BreakPoint[i])==curtype && curtype!=3){
                output<<'\t'<<get<1>(BreakPoint[i])<<'\t'<<curid<<'\t'<<SimpleSVTypeString[get<3>(BreakPoint[i])]<<endl;
            }
            else if (curtype==3){
                i--;
                output<<"\t.\t"<<curid<<'\t'<<"TRA"<<endl;
            }
            i++;
        } while(i<sizeBP);
        output.flush();
        output.close();
    }
    else
        cout<<"Cannot open output file \""<<outputfile<<"\"\n";
};

void SV_t::WritenewSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile){
    if(!updated)
        Update2NewPos(RefLength);
    ofstream output(outputfile, ios::out);
    output<<"## Type\tChr1\tPos1\tDT1\tChr2\tPos2\tDT2\tID"<<endl;
    for(vector<SimpleSV_t>::iterator it=vSimpleSV.begin(); it!=vSimpleSV.end(); it++){
        output<<SimpleSVTypeString[it->Type]<<'\t'<<RefName[it->RefID]<<'\t'<<it->StartPos<<"\t.\t.\t"<<it->EndPos<<"\t.\t"<<it->ID<<"\n";
    }
    for(vector<TRA_t>::iterator it=vTRA.begin(); it!=vTRA.end(); it++){
        output<<"TRA\t"<<RefName[it->Ref1]<<'\t'<<it->Pos1<<'\t'<<it->DT1<<'\t'<<RefName[it->Ref2]<<'\t'<<it->Pos2<<'\t'<<it->DT2<<'\t'<<it->ID<<'\n';
    }
    output.close();
};

void SV_t::WriteoldSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile){
    if(updated)
        Update2NewPos(RefLength);
    ofstream output(outputfile, ios::out);
    output<<"## Type\tChr1\tPos1\tDT1\tChr2\tPos2\tDT2\tID"<<endl;
    for(vector<SimpleSV_t>::iterator it=vSimpleSV.begin(); it!=vSimpleSV.end(); it++){
        output<<SimpleSVTypeString[it->Type]<<'\t'<<RefName[it->RefID]<<'\t'<<it->StartPos<<"\t.\t.\t"<<it->EndPos<<"\t.\t"<<it->ID<<"\n";
    }
    for(vector<TRA_t>::iterator it=vTRA.begin(); it!=vTRA.end(); it++){
        output<<"TRA\t"<<RefName[it->Ref1]<<'\t'<<it->Pos1<<'\t'<<it->DT1<<'\t'<<RefName[it->Ref2]<<'\t'<<it->Pos2<<'\t'<<it->DT2<<'\t'<<it->ID<<'\n';
    }
    output.close();
};

void SV_t::WriteFilterednewSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile){
    if(vSimpleReads.size()==0 || vTRAReads.size()==0){
        cout<<"Error: it is not calculated yet whether SVs are covered by reads!"<<endl;
        return;
    }
    if(!updated)
        Update2NewPos(RefLength);
    ofstream output(outputfile, ios::out);
    output<<"## Type\tChr1\tPos1\tDT1\tChr2\tPos2\tDT2\tID"<<endl;
    int count=0;
    for(vector<SimpleSV_t>::iterator it=vSimpleSV.begin(); it!=vSimpleSV.end(); it++){
        if(vSimpleReads[count]){
            bool flag=true;
            // filter out single duplications and deletions, but keep insertions (transposons / transpositions)
            if(vSimpleSV[count].Type==SimpleSVType::INS && ((count==0)?true:(vSimpleSV[count-1].ID!=vSimpleSV[count].ID)) && ((count==vSimpleSV.size()-1)?true:(vSimpleSV[count+1].ID!=vSimpleSV[count].ID)))
                flag=false;
            else if(vSimpleSV[count].Type==SimpleSVType::DEL && ((count==0)?true:(vSimpleSV[count-1].ID!=vSimpleSV[count].ID)) && ((count==vSimpleSV.size()-1)?true:(vSimpleSV[count+1].ID!=vSimpleSV[count].ID)))
                flag=false;
            if(flag)
                output<<SimpleSVTypeString[it->Type]<<'\t'<<RefName[it->RefID]<<'\t'<<it->StartPos<<"\t.\t.\t"<<it->EndPos<<"\t.\t"<<it->ID<<"\n";
        }
        count++;
    }
    count=0;
    for(vector<TRA_t>::iterator it=vTRA.begin(); it!=vTRA.end(); it++){
        if(vTRAReads[count++])
            output<<"TRA\t"<<RefName[it->Ref1]<<'\t'<<it->Pos1<<'\t'<<it->DT1<<'\t'<<RefName[it->Ref2]<<'\t'<<it->Pos2<<'\t'<<it->DT2<<'\t'<<it->ID<<'\n';
    }
    output.close();
};

void SV_t::WriteFilteredoldSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile){
    if(vSimpleReads.size()==0 || vTRAReads.size()==0){
        cout<<"Error: it is not calculated yet whether SVs are covered by reads!"<<endl;
        return;
    }
    if(updated)
        Update2NewPos(RefLength);
    ofstream output(outputfile, ios::out);
    output<<"## Type\tChr1\tPos1\tDT1\tChr2\tPos2\tDT2\tID"<<endl;
    int count=0;
    for(vector<SimpleSV_t>::iterator it=vSimpleSV.begin(); it!=vSimpleSV.end(); it++){
        if(vSimpleReads[count]){
            bool flag=true;
            // filter out single duplications and deletions, but keep insertions (transposons / transpositions)
            if(vSimpleSV[count].Type==SimpleSVType::INS && ((count==0)?true:(vSimpleSV[count-1].ID!=vSimpleSV[count].ID)) && ((count==vSimpleSV.size()-1)?true:(vSimpleSV[count+1].ID!=vSimpleSV[count].ID)))
                flag=false;
            else if(vSimpleSV[count].Type==SimpleSVType::DEL && ((count==0)?true:(vSimpleSV[count-1].ID!=vSimpleSV[count].ID)) && ((count==vSimpleSV.size()-1)?true:(vSimpleSV[count+1].ID!=vSimpleSV[count].ID)))
                flag=false;
            if(flag)
                output<<SimpleSVTypeString[it->Type]<<'\t'<<RefName[it->RefID]<<'\t'<<it->StartPos<<"\t.\t.\t"<<it->EndPos<<"\t.\t"<<it->ID<<"\n";
        }
        count++;
    }
    count=0;
    for(vector<TRA_t>::iterator it=vTRA.begin(); it!=vTRA.end(); it++){
        if(vTRAReads[count++])
            output<<"TRA\t"<<RefName[it->Ref1]<<'\t'<<it->Pos1<<'\t'<<it->DT1<<'\t'<<RefName[it->Ref2]<<'\t'<<it->Pos2<<'\t'<<it->DT2<<'\t'<<it->ID<<'\n';
    }
    output.close();
};

bool SV_t::IsIntersect(SV_t svrhs){
    bool Intersect=false;
    sort(vSimpleSV.begin(), vSimpleSV.end());
    sort(svrhs.vSimpleSV.begin(), svrhs.vSimpleSV.end());
    vector<SimpleSV_t>::iterator it1=vSimpleSV.begin(), it2=svrhs.vSimpleSV.begin();
    for(;it1!=vSimpleSV.end() && it2!=svrhs.vSimpleSV.end();){
        if((*it1)<(*it2)){
            if(it1->RefID==it2->RefID && it1->EndPos>it2->StartPos){
                Intersect=true;
                cout<<'\t'<<it1->Print();
                cout<<'\t'<<it2->Print();
                break;
            }
            it1++;
        }
        else{
            if(it1->RefID==it2->RefID && it1->StartPos<it2->EndPos){
                Intersect=true;
                cout<<'\t'<<it1->Print();
                cout<<'\t'<<it2->Print();
                break;
            }
            it2++;
        }
    }
    return Intersect;
};

vector< pair<int,int> > Convert2GenomeCoord(pair<int,int> interval, Transcript_t& currenttrans)
{
    vector< pair<int,int> > newintervals;
    if(currenttrans.Strand=='+'){
        int newintvlstart=interval.first;
        int exonidx=0;
        int totallen=0;
        while(newintvlstart<interval.second){
            while(totallen+currenttrans.vExon[exonidx].EndPos-currenttrans.vExon[exonidx].StartPos<=newintvlstart){
                totallen+=currenttrans.vExon[exonidx].EndPos-currenttrans.vExon[exonidx].StartPos;
                exonidx++;
            }
            assert(exonidx < (int)currenttrans.vExon.size());
            int step=min(interval.second-newintvlstart, currenttrans.vExon[exonidx].EndPos-currenttrans.vExon[exonidx].StartPos-(newintvlstart-totallen));
            int newend=newintvlstart+step;
            newintervals.push_back(make_pair(currenttrans.vExon[exonidx].StartPos+newintvlstart-totallen, currenttrans.vExon[exonidx].StartPos+newend-totallen));
            assert(newintervals.back().first<=currenttrans.vExon[exonidx].EndPos && newintervals.back().second<=currenttrans.vExon[exonidx].EndPos);
            newintvlstart=newend;
        }
    }
    else{
        int newintvlstart=interval.first;
        int exonidx=(int)currenttrans.vExon.size()-1;
        int totallen=0;
        while(newintvlstart<interval.second){
            while(totallen+currenttrans.vExon[exonidx].EndPos-currenttrans.vExon[exonidx].StartPos<=newintvlstart){
                totallen+=currenttrans.vExon[exonidx].EndPos-currenttrans.vExon[exonidx].StartPos;
                exonidx--;
            }
            if(exonidx<0)
                cout<<"watch here\n";
            assert(exonidx>=0);
            int step=min(interval.second-newintvlstart, currenttrans.vExon[exonidx].EndPos-currenttrans.vExon[exonidx].StartPos-(newintvlstart-totallen));
            int newend=newintvlstart+step;
            newintervals.push_back(make_pair(currenttrans.vExon[exonidx].EndPos-(newintvlstart-totallen), currenttrans.vExon[exonidx].EndPos-(newend-totallen)));
            assert(newintervals.back().first>=currenttrans.vExon[exonidx].StartPos && newintervals.back().second>=currenttrans.vExon[exonidx].StartPos);
            newintvlstart=newend;
        }
        reverse(newintervals.begin(), newintervals.end());
    }
    return newintervals;
};

void SV_t::IsCoveredbyPolyester(string fastafile, vector<Transcript_t>& vTrans, map<int,int>& RefLength, vector<string>& TransName, bool is_old_coordicate, int threshold){
    if(is_old_coordicate && updated)
        Update2NewPos(RefLength);
    else if(!is_old_coordicate && !updated)
        Update2NewPos(RefLength);
    vSimpleReads.clear(); vSimpleReads.resize(vSimpleSV.size(), false);
    vTRAReads.clear(); vTRAReads.resize(vTRA.size(), false);
    vector<int> countSimpleReads(vSimpleSV.size(), 0);
    vector<int> countTRAReads(vTRA.size(), 0);
    // learn read length from simulated fasta file
    int readlen=0;

    // process transcripts vector for quick access with transcript name
    map<string,int> TransIndexMap;
    for(int i=0; i<vTrans.size(); i++)
        TransIndexMap[vTrans[i].TransID]=i;

    ifstream input(fastafile);
    string line;
    while(getline(input,line)){
        if(line[0]!='>')
            continue;
        vector<string> strs;
        boost::split(strs, line, boost::is_any_of(";"));

        size_t s=strs[0].find_first_of('/');
        string transname=strs[0].substr(s+1);
        Transcript_t& currenttrans=vTrans[TransIndexMap[transname]];
        int readchr=currenttrans.Chr;
        // get mate 1 and 2 intervals in transcript coordinate
        pair<int,int> mate1_interval_ori;
        pair<int,int> mate2_interval_ori;
        assert(strs[1].substr(0, 5)=="mate1" && strs[2].substr(0,5)=="mate2");
        if(strs[1][5]==':'){
            strs[1]=strs[1].substr(6);
            vector<string> mate1str;
            boost::split(mate1str, strs[1], boost::is_any_of("-"));
            assert(mate1str.size()>1);
            mate1_interval_ori=make_pair(stoi(mate1str[0])-1, stoi(mate1str[1]));
            if(readlen==0)
                readlen = mate1_interval_ori.second-mate1_interval_ori.first;
        }
        else if(strs[1].substr(5,6)=="Start:"){
            strs[1]=strs[1].substr(11);
            mate1_interval_ori=make_pair(stoi(strs[1])-1, stoi(strs[1])-1+readlen);
            // check corresponding transcript length
            int translen=0;
            for(int i=0; i<currenttrans.vExon.size(); i++)
                translen+=currenttrans.vExon[i].EndPos-currenttrans.vExon[i].StartPos;
            if(translen<readlen)
                mate1_interval_ori.second=stoi(strs[1])-1+translen;
        }
        else
            cout<<"Invalid format"<<line<<endl;
        if(strs[2][5]==':'){
            strs[2]=strs[2].substr(6);
            vector<string> mate2str;
            boost::split(mate2str, strs[2], boost::is_any_of("-"));
            assert(mate2str.size()>1);
            mate2_interval_ori=make_pair(stoi(mate2str[0])-1, stoi(mate2str[1]));
            if(readlen==0)
                readlen = mate2_interval_ori.second-mate2_interval_ori.first;
        }
        else if(strs[2].substr(5,6)=="Start:"){
            strs[2]=strs[2].substr(11);
            mate2_interval_ori=make_pair(stoi(strs[2])-1, stoi(strs[2])-1+readlen);
            // check corresponding transcript length
            int translen=0;
            for(int i=0; i<currenttrans.vExon.size(); i++)
                translen+=currenttrans.vExon[i].EndPos-currenttrans.vExon[i].StartPos;
            if(translen<readlen)
                mate2_interval_ori.second=stoi(strs[1])-1+translen;
        }
        else
            cout<<"Invalid format"<<line<<endl;

        //  convert mate 1 and 2 positions into genomic coordinate
        vector< pair<int,int> > mate1_interval_conv=Convert2GenomeCoord(mate1_interval_ori, currenttrans);
        vector< pair<int,int> > mate2_interval_conv=Convert2GenomeCoord(mate2_interval_ori, currenttrans);
        int numfirstread=mate1_interval_conv.size();
        vector<int> Positions;
        vector<int> Lengths;
        for(vector< pair<int,int> >::iterator it=mate1_interval_conv.begin(); it!=mate1_interval_conv.end(); it++){
            Positions.push_back(it->first);
            Lengths.push_back(it->second - it->first);
        }
        for(vector< pair<int,int> >::iterator it=mate2_interval_conv.begin(); it!=mate2_interval_conv.end(); it++){
            if(currenttrans.Strand=='+'){
                Positions.push_back(it->first);
                Lengths.push_back(it->second - it->first);
            }
            else{
                Positions.insert(Positions.begin(), it->first);
                Lengths.insert(Lengths.begin(), it->second - it->first);
            }
        }

        // variables to keep track of which TSVs are supported
        vector<int> SimpleIdxList;
        vector<int> TRAIdxList;
        bool SupportingBP=false;

        // check simple SVs
        for(int i=0; i<vSimpleSV.size(); i++){
            if(strs[0].substr(0,12)==">read2307089" && vSimpleSV[i].ID==327)
                cout<<"watch here\n";
            if(readchr!=vSimpleSV[i].RefID)
                continue;
            for(int j=0; j<Positions.size(); j++)
                if((Positions[j]<vSimpleSV[i].StartPos-10 && Positions[j]+Lengths[j]>vSimpleSV[i].StartPos+10) || (Positions[j]<vSimpleSV[i].EndPos-10 && Positions[j]+Lengths[j]>vSimpleSV[i].EndPos+10)){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            for(int j=0; j<Positions.size()-1; j++){
                if(Positions[j]+Lengths[j]<vSimpleSV[i].StartPos+10 && Positions[j+1]>vSimpleSV[i].StartPos-10 && Positions[j+1]+Lengths[j+1]<vSimpleSV[i].EndPos+10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
                else if(Positions[j]>vSimpleSV[i].StartPos-10 && Positions[j]+Lengths[j]<vSimpleSV[i].EndPos+10 && Positions[j+1]>vSimpleSV[i].EndPos-10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            }
            if(numfirstread<Positions.size() && Positions[0]<Positions[numfirstread]){
                if(Positions[0]+Lengths[0]<vSimpleSV[i].StartPos+10 && Positions.back()>vSimpleSV[i].StartPos-10 && Positions.back()+Lengths.back()<vSimpleSV[i].EndPos+10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
                else if(Positions[0]>vSimpleSV[i].StartPos-10 && Positions[0]+Lengths[0]<vSimpleSV[i].EndPos+10 && Positions.back()>vSimpleSV[i].EndPos-10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            }
            else if(numfirstread<Positions.size()){
                if(Positions[numfirstread]+Lengths[numfirstread]<vSimpleSV[i].StartPos+10 && Positions[numfirstread-1]>vSimpleSV[i].StartPos-10 && Positions[numfirstread-1]+Lengths[numfirstread-1]<vSimpleSV[i].EndPos+10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
                else if(Positions[numfirstread]>vSimpleSV[i].StartPos-10 && Positions[numfirstread]+Lengths[numfirstread]<vSimpleSV[i].EndPos+10 && Positions[numfirstread-1]>vSimpleSV[i].EndPos-10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            }
        }
        // check TRAs
        for(int i=0; i<vTRA.size(); i++){
            if(readchr!=vTRA[i].Ref1 && readchr!=vTRA[i].Ref2)
                continue;
            for(int j=0; j<Positions.size(); j++)
                if((readchr==vTRA[i].Ref1 && Positions[j]<vTRA[i].Pos1-10 && Positions[j]+Lengths[j]>vTRA[i].Pos1+10) || (readchr==vTRA[i].Ref2 && Positions[j]<vTRA[i].Pos2-10 && Positions[j]+Lengths[j]>vTRA[i].Pos2+10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            for(int j=0; j<Positions.size()-1; j++)
                if((readchr==vTRA[i].Ref1 && Positions[j]+Lengths[j]<vTRA[i].Pos1+10 && Positions[j+1]>vTRA[i].Pos1-10) || (readchr==vTRA[i].Ref2 && Positions[j]+Lengths[j]<vTRA[i].Pos2+10 && Positions[j+1]>vTRA[i].Pos2-10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            if(numfirstread<Positions.size() && Positions[0]<Positions[numfirstread]){
                if((readchr==vTRA[i].Ref1 && Positions[0]+Lengths[0]<vTRA[i].Pos1+10 && Positions.back()>vTRA[i].Pos1-10) || (readchr==vTRA[i].Ref2 && Positions[0]+Lengths[0]<vTRA[i].Pos2+10 && Positions.back()>vTRA[i].Pos2-10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            }
            else if(numfirstread<Positions.size()){
                if((readchr==vTRA[i].Ref1 && Positions[numfirstread]+Lengths[numfirstread]<vTRA[i].Pos1+10 && Positions[numfirstread-1]>vTRA[i].Pos1-10) || (readchr==vTRA[i].Ref2 && Positions[numfirstread]+Lengths[numfirstread]<vTRA[i].Pos2+10 && Positions[numfirstread-1]>vTRA[i].Pos2-10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            }
        }
        if(SupportingBP)
            TransName.push_back(transname);
        // add to count
        sort(SimpleIdxList.begin(), SimpleIdxList.end());
        SimpleIdxList.resize(distance(SimpleIdxList.begin(), unique(SimpleIdxList.begin(), SimpleIdxList.end())));
        sort(TRAIdxList.begin(), TRAIdxList.end());
        TRAIdxList.resize(distance(TRAIdxList.begin(), unique(TRAIdxList.begin(), TRAIdxList.end())));
        for(int i=0; i<SimpleIdxList.size(); i++){
            countSimpleReads[SimpleIdxList[i]]++;
        }
        for(int i=0; i<TRAIdxList.size(); i++){
            countTRAReads[TRAIdxList[i]]++;
        }
    }
    input.close();

    // mark TSVs to be supported
    for(int i=0; i<countSimpleReads.size(); i++)
        if(countSimpleReads[i]>threshold)
            vSimpleReads[i]=true;
    for(int i=0; i<countTRAReads.size(); i++)
        if(countTRAReads[i]>threshold)
            vTRAReads[i]=true;
    // for same SV ID, if one is covered by reads, mark the others covered.
    int idcount=0;
    for(int i=0; i<vSimpleSV.size(); i++)
        if(idcount<vSimpleSV[i].ID)
            idcount=vSimpleSV[i].ID;
    for(int i=0; i<vTRA.size(); i++)
        if(idcount<vTRA[i].ID)
            idcount=vTRA[i].ID;
    for(int i=0; i<idcount; i++){
        bool iscovered=false;
        for(int j=0; j<vSimpleSV.size(); j++)
            if(vSimpleSV[j].ID==i && vSimpleReads[j])
                iscovered=true;
        for(int j=0; j<vTRA.size(); j++)
            if(vTRA[j].ID==i && vTRAReads[j])
                iscovered=true;
        if(iscovered){
            for(int j=0; j<vSimpleSV.size(); j++)
                if(vSimpleSV[j].ID==i)
                    vSimpleReads[j]=true;
            for(int j=0; j<vTRA.size(); j++)
                if(vTRA[j].ID==i)
                    vTRAReads[j]=true;
        }
    }
    // trim transname
    sort(TransName.begin(), TransName.end());
    vector<string>::iterator endittrans=unique(TransName.begin(), TransName.end());
    TransName.resize(distance(TransName.begin(), endittrans));
};

void SV_t::IsCoveredbyReads(string bedfile, map<string,int>& RefTable, map<int,int>& RefLength, vector<string>& TransName, bool is_old_coordicate, int threshold){
    if(is_old_coordicate && updated)
        Update2NewPos(RefLength);
    else if(!is_old_coordicate && !updated)
        Update2NewPos(RefLength);
    vSimpleReads.clear(); vSimpleReads.resize(vSimpleSV.size(), false);
    vTRAReads.clear(); vTRAReads.resize(vTRA.size(), false);
    vector< vector< vector<int> > > countSimpleReads(vSimpleSV.size());
    vector< vector< vector<int> > > countTRAReads(vTRA.size());
    vector< vector<string> > readnameSimpleReads(vSimpleSV.size());
    vector< vector<string> > readnameTRAReads(vTRA.size());

    int actualsupport=0;

    clock_t starttime=clock();
    double duration;
    ifstream input(bedfile);
    string line;
    bool can_getline=getline(input, line);
    int readnum=0;
    while(can_getline){
        vector<int> Positions, Lengths;
        int numfirstread=0, readchr=0;
        vector<string> read1str, read2str, blocklen, offset;
        boost::split(read1str, line, boost::is_any_of("\t"));
        can_getline=getline(input, line);
        boost::split(read2str, line, boost::is_any_of("\t"));
        can_getline=getline(input, line);
        readnum++;
        if(readnum%1000000==0){
            duration=1.0*(clock()-starttime)/CLOCKS_PER_SEC;
            cout<<readnum<<"\tused time "<<duration<<endl;
        }
        readchr=RefTable[read1str[0]];
        // read 1
        if(read1str.size()>10){
            if(read1str[0]!="polyA"){
                boost::split(blocklen, read1str[10], boost::is_any_of(","));
                boost::split(offset, read1str[11], boost::is_any_of(","));
                for(int i=0; i<offset.size(); i++){
                    Positions.push_back(stoi(read1str[1])+stoi(offset[i]));
                    Lengths.push_back(stoi(blocklen[i]));
                }
            }
            numfirstread=blocklen.size();
        }
        else{
            Positions.push_back(stoi(read1str[1]));
            Lengths.push_back(stoi(read1str[2])-stoi(read1str[1]));
            numfirstread=1;
        }
        // read 2
        if(read2str.size()>10){
            if(read2str[0]!="polyA"){
                boost::split(blocklen, read2str[10], boost::is_any_of(","));
                boost::split(offset, read2str[11], boost::is_any_of(","));
                for(int i=0; i<offset.size(); i++){
                    Positions.push_back(stoi(read2str[1])+stoi(offset[i]));
                    Lengths.push_back(stoi(blocklen[i]));
                }
            }
        }
        else{
            Positions.push_back(stoi(read2str[1]));
            Lengths.push_back(stoi(read2str[2])-stoi(read2str[1]));
        }
        // skip reads where both ends are polyA
        if(Positions.size()==0)
            continue;
        // skip reads that have too many segments
        if(numfirstread>2)
            continue;
        vector<int> SimpleIdxList;
        vector<int> TRAIdxList;
        // find in SV lists
        bool SupportingBP=false;
        for(int i=0; i<vSimpleSV.size(); i++){
            if(readchr!=vSimpleSV[i].RefID)
                continue;
            for(int j=0; j<Positions.size(); j++)
                if((Positions[j]<vSimpleSV[i].StartPos-10 && Positions[j]+Lengths[j]>vSimpleSV[i].StartPos+10) || (Positions[j]<vSimpleSV[i].EndPos-10 && Positions[j]+Lengths[j]>vSimpleSV[i].EndPos+10)){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            for(int j=0; j<Positions.size()-1; j++){
                if(Positions[j]+Lengths[j]<vSimpleSV[i].StartPos+10 && Positions[j+1]>vSimpleSV[i].StartPos-10 && Positions[j+1]+Lengths[j+1]<vSimpleSV[i].EndPos+10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
                else if(Positions[j]>vSimpleSV[i].StartPos-10 && Positions[j]+Lengths[j]<vSimpleSV[i].EndPos+10 && Positions[j+1]>vSimpleSV[i].EndPos-10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            }
            if(numfirstread<Positions.size() && Positions[0]<Positions[numfirstread]){
                if(Positions[0]+Lengths[0]<vSimpleSV[i].StartPos+10 && Positions.back()>vSimpleSV[i].StartPos-10 && Positions.back()+Lengths.back()<vSimpleSV[i].EndPos+10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
                else if(Positions[0]>vSimpleSV[i].StartPos-10 && Positions[0]+Lengths[0]<vSimpleSV[i].EndPos+10 && Positions.back()>vSimpleSV[i].EndPos-10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            }
            else if(numfirstread<Positions.size()){
                if(Positions[numfirstread]+Lengths[numfirstread]<vSimpleSV[i].StartPos+10 && Positions[numfirstread-1]>vSimpleSV[i].StartPos-10 && Positions[numfirstread-1]+Lengths[numfirstread-1]<vSimpleSV[i].EndPos+10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
                else if(Positions[numfirstread]>vSimpleSV[i].StartPos-10 && Positions[numfirstread]+Lengths[numfirstread]<vSimpleSV[i].EndPos+10 && Positions[numfirstread-1]>vSimpleSV[i].EndPos-10){
                    SimpleIdxList.push_back(i);
                    assert(i<vSimpleSV.size());
                    SupportingBP=true;
                }
            }
        }
        for(int i=0; i<vTRA.size(); i++){
            if(readchr!=vTRA[i].Ref1 && readchr!=vTRA[i].Ref2)
                continue;
            for(int j=0; j<Positions.size(); j++)
                if((readchr==vTRA[i].Ref1 && Positions[j]<vTRA[i].Pos1-10 && Positions[j]+Lengths[j]>vTRA[i].Pos1+10) || (readchr==vTRA[i].Ref2 && Positions[j]<vTRA[i].Pos2-10 && Positions[j]+Lengths[j]>vTRA[i].Pos2+10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            for(int j=0; j<Positions.size()-1; j++)
                if((readchr==vTRA[i].Ref1 && Positions[j]+Lengths[j]<vTRA[i].Pos1+10 && Positions[j+1]>vTRA[i].Pos1-10) || (readchr==vTRA[i].Ref2 && Positions[j]+Lengths[j]<vTRA[i].Pos2+10 && Positions[j+1]>vTRA[i].Pos2-10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            if(numfirstread<Positions.size() && Positions[0]<Positions[numfirstread]){
                if((readchr==vTRA[i].Ref1 && Positions[0]+Lengths[0]<vTRA[i].Pos1+10 && Positions.back()>vTRA[i].Pos1-10) || (readchr==vTRA[i].Ref2 && Positions[0]+Lengths[0]<vTRA[i].Pos2+10 && Positions.back()>vTRA[i].Pos2-10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            }
            else if(numfirstread<Positions.size()){
                if((readchr==vTRA[i].Ref1 && Positions[numfirstread]+Lengths[numfirstread]<vTRA[i].Pos1+10 && Positions[numfirstread-1]>vTRA[i].Pos1-10) || (readchr==vTRA[i].Ref2 && Positions[numfirstread]+Lengths[numfirstread]<vTRA[i].Pos2+10 && Positions[numfirstread-1]>vTRA[i].Pos2-10)){
                    TRAIdxList.push_back(i);
                    assert(i<vTRA.size());
                    SupportingBP=true;
                }
            }
        }
        if(SupportingBP){
            vector<string> readnameinfo;
            boost::split(readnameinfo, read1str[3], boost::is_any_of(":"));
            TransName.push_back(readnameinfo[2]);
        }
        vector<int>::iterator endit;
        sort(SimpleIdxList.begin(), SimpleIdxList.end());
        endit=unique(SimpleIdxList.begin(), SimpleIdxList.end());
        SimpleIdxList.resize(distance(SimpleIdxList.begin(), endit));
        sort(TRAIdxList.begin(), TRAIdxList.end());
        endit=unique(TRAIdxList.begin(), TRAIdxList.end());
        TRAIdxList.resize(distance(TRAIdxList.begin(), endit));
        for(int i=0; i<SimpleIdxList.size(); i++){
            countSimpleReads[SimpleIdxList[i]].push_back(Positions);
            readnameSimpleReads[SimpleIdxList[i]].push_back(read1str[3].substr(0, (int)read1str[3].size()-2));
        }
        for(int i=0; i<TRAIdxList.size(); i++){
            countTRAReads[TRAIdxList[i]].push_back(Positions);
            readnameTRAReads[TRAIdxList[i]].push_back(read1str[3].substr(0, (int)read1str[3].size()-2));
        }
    }
    input.close();

    sort(TransName.begin(), TransName.end());
    vector<string>::iterator endittrans=unique(TransName.begin(), TransName.end());
    TransName.resize(distance(TransName.begin(), endittrans));

    vector<string> ValidReadName;
    // for same SV ID, if one is covered by reads, mark the others covered.
    for(int i=0; i<countSimpleReads.size(); i++){
        int countsupport=0;
        for(int j=0; j<countSimpleReads[i].size(); j++){
            // test this read position vector is same as previous ones
            bool duplicate=false;
            for(int k=j-1; k>=0; k--){
                bool same=true;
                if(countSimpleReads[i][j].size()!=countSimpleReads[i][k].size())
                    same=false;
                else{
                    for(int s=0; s<countSimpleReads[i][j].size(); s++)
                        if(countSimpleReads[i][j][s]!=countSimpleReads[i][k][s])
                            same=false;
                }
                if(same){
                    duplicate=true;
                    break;
                }
                if(duplicate)
                    break;
            }
            if(!duplicate){
                countsupport++;
                ValidReadName.push_back(readnameSimpleReads[i][j]);
            }
        }
        if(countsupport>threshold)
            vSimpleReads[i]=true;
        actualsupport+=countsupport;
    }
    for(int i=0; i<countTRAReads.size(); i++){
        int countsupport=0;
        for(int j=0; j<countTRAReads[i].size(); j++){
            // test this read position vector is same as previous ones
            bool duplicate=false;
            for(int k=j-1; k>=0; k--){
                bool same=true;
                if(countTRAReads[i][j].size()!=countTRAReads[i][k].size())
                    same=false;
                else{
                    for(int s=0; s<countTRAReads[i][j].size(); s++)
                        if(countTRAReads[i][j][s]!=countTRAReads[i][k][s])
                            same=false;
                }
                if(same){
                    duplicate=true;
                    break;
                }
                if(duplicate)
                    break;
            }
            if(!duplicate){
                countsupport++;
                ValidReadName.push_back(readnameTRAReads[i][j]);
            }
        }
        if(countsupport>threshold)
            vTRAReads[i]=true;
        actualsupport+=countsupport;
    }
    int idcount=0;
    for(int i=0; i<vSimpleSV.size(); i++)
        if(idcount<vSimpleSV[i].ID)
            idcount=vSimpleSV[i].ID;
    for(int i=0; i<vTRA.size(); i++)
        if(idcount<vTRA[i].ID)
            idcount=vTRA[i].ID;
    for(int i=0; i<idcount; i++){
        bool iscovered=false;
        for(int j=0; j<vSimpleSV.size(); j++)
            if(vSimpleSV[j].ID==i && vSimpleReads[j])
                iscovered=true;
        for(int j=0; j<vTRA.size(); j++)
            if(vTRA[j].ID==i && vTRAReads[j])
                iscovered=true;
        if(iscovered){
            for(int j=0; j<vSimpleSV.size(); j++)
                if(vSimpleSV[j].ID==i)
                    vSimpleReads[j]=true;
            for(int j=0; j<vTRA.size(); j++)
                if(vTRA[j].ID==i)
                    vTRAReads[j]=true;
        }
    }

    cout<<"Actual number of TSV supporting reads = "<<actualsupport<<endl;
    // sort(ValidReadName.begin(), ValidReadName.end());
    // ValidReadName.resize(distance(ValidReadName.begin(), unique(ValidReadName.begin(), ValidReadName.end())));
    // for(int i=0; i<ValidReadName.size(); i++)
    //     cout<<ValidReadName[i]<<endl;
};
