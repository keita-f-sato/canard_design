#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
#define g 9.80665
#define lo 1.253
#define pi 3.14159265358979323846
#define e 0.8
#define nyu 0.0000146

double CLw,CLt,ARw,ARt,m,x,v,as,Rew,Ret,CDw,CDt,CDwt,CL2w,CL2t;
double mkg,mN,Lwkg,Ltkg,Lwg,Ltg,a0,alphaw,alphat,daw,dat,LNw,LNt;
double mas,cw,ct,cwm,ctm,lw,lt,bw,bt,Sw,St,Lt,Lw,ldt,arw,art;
double aw,at,Vh,hn,hnw,h,Cmalpha,CLalpha,dalphade,a0w,a0t,hd,lt1;
int i,i2,i3,NAw,NAt;

class Airplane{
public:
	  int iuputdata();
	  int analysisdata();
	  int outputdata();
};

int Airplane::iuputdata(){
    cout<<"-------���œK���݌v-------"<<endl;
    cout<<"-------�嗃-------"<<endl;
    cout<<"�嗃�̊e�f�[�^����͂��Ă�������"<<endl;
    cout<<"NACA";
    cin>>NAw;
    cout<<"�嗃�A�X�y�N�g��:";
    cin>>ARw;
    cout<<"�嗃���g�͊p[deg]:";
    cin>>a0w;
    cout<<"�嗃�}�p(��)[deg]";
    cin>>alphaw;
    cout<<endl;

    cout<<"-------�����-------"<<endl;
    cout<<"������̊e�f�[�^����͂��Ă�������"<<endl;
    cout<<"NACA";
    cin>>NAt;

    cout<<"������A�X�y�N�g��:";
    cin>>ARt;
    cout<<"��������g�͊p[deg]:";
    cin>>a0t;
    cout<<"�����}�p(��)[deg]";
    cin>>alphat;
    cout<<endl;

    cout<<"-------����-------"<<endl;
    cout<<"�@�̏�������͂��Ă�������"<<endl;
    cout<<"�@�̏d��[g]�F";
    cin>>m;
    cout<<"���q���x[m/s]�F";
    cin>>v;
    cout<<"-------�݌v��-------"<<endl;
    cout<<"���S���F";
    cin>>as;
    cout<<"�O��g�͔z��Lt�FLw=1�F";
    cin>>x;
    return(0);

}

int Airplane::analysisdata(){
    mkg=m/1000;
    mN=mkg*g;
    mas=mN*as;
    Lt=mas/(1+x);
    Lw=x*Lt;

    Lwkg=Lw/g;
    Ltkg=Lt/g;
    Lwg=Lwkg*1000;
    Ltg=Ltkg*1000;

    a0=2*pi;
    arw=alphaw*(180/pi);
    art=alphat*(180/pi);
    CL2w=a0*arw;
    CL2t=a0*art;

    daw=CL2w/(pi*ARw);
    dat=CL2t/(pi*ARt);
    aw=daw/(1+(daw/(pi*ARw*e)));
    at=dat/(1+(dat/(pi*ARt*e)));

    CLw=daw*(-1*a0w)/(1+(daw/pi*ARw*e));
    CLt=dat*(-1*a0t)/(1+(daw/pi*ARw*e));
    cwm=sqrt((2*x*mas)/(lo*v*v*CLw*ARw*(1+x)));
    ctm=sqrt((2*mas)/(lo*v*v*CLt*ARt*(1+x)));
    cw=cwm*1000;
    ct=ctm*1000;
    bw=cw*ARw;
    bt=ct*ARt;
    Sw=cw*bw;
    St=ct*bt;



    dalphade=(2*aw)/(pi*ARw);

    CLalpha=aw*(1+(at/aw)*(St/Sw)*(1-dalphade));
    hnw=0.25;
    lt=0.5*((Lw*cw*hnw)/Lt)+(aw/at)*(1/(1-dalphade))*(Sw*cw/St*hnw);
    hn=hnw+(1/CLalpha)*(at*(1-dalphade)*lt*(St/(Sw*ct)));
    lw=(hn/2)*cw+0.25*cw;
    lt1=Lw*lw/Lt;

    Rew=v*cwm/nyu;
    Ret=v*ctm/nyu;

    cout<<"-------��́��݌v����-------"<<endl;
    cout<<"-------���ʐ�-------"<<endl;
    cout<<"�嗃�S��"<<bw<<endl;
    cout<<"�嗃�R�[�h��"<<cw<<endl;
    cout<<"�嗃�ʐ�"<<Sw<<endl;
    cout<<"�����S��"<<bt<<endl;
    cout<<"�����R�[�h��"<<ct<<endl;
    cout<<"�����ʐ�"<<St<<endl;
    cout<<"-------�����g��-------"<<endl;
    cout<<"�嗃���C�m���Y��[-]"<<Rew<<endl;
    cout<<"�嗃�g�͌W��[-]"<<CLw<<endl;
    cout<<"�嗃�����g��[g]"<<Lwg<<endl;
    cout<<"�嗃�g�͌W��[-]"<<CLw<<endl;
    cout<<"�������C�m���Y��[-]"<<Ret<<endl;
    cout<<"�����g�͌W��[-]"<<CLt<<endl;
    cout<<"���������g��[g]"<<Ltg<<endl;
    cout<<"���v�����g��[g]"<<Lwg+Ltg<<endl;
    cout<<"-------����݌v-------"<<endl;
    cout<<"CL��:"<<CLalpha<<endl;
    cout<<"�݌v�d�S�ʒu[%]"<<hn/2<<endl;
    cout<<"�嗃��͒��S����̏d�S�ʒu[mm]"<<lw<<endl;
    cout<<"������͒��S����̏d�S�ʒu[mm]"<<lt1<<endl;
    cout<<"�嗃�O������̏d�S�ʒu[mm]"<<cw*(hn/2)<<endl;
    cout<<"�����O������̏d�S�ʒu[mm]"<<lt1+ct*hnw<<endl;

   return(0);
}

int Airplane::outputdata(){

    ofstream of;
    of.open("Airplane.txt",ios::app);
    of<<endl;
    of<<"-------���œK���݌v-------"<<endl;
    of<<"-------�嗃-------"<<endl;
    of<<"NACA"<<NAw<<endl;
    of<<"�嗃�A�X�y�N�g��:"<<ARw<<endl;
    of<<"�嗃���g�͊p"<<a0w<<"[deg]"<<endl;
    of<<"�嗃�}�p(��)"<<alphaw<<"[deg]"<<endl;
    of<<endl;
    of<<"-------�����-------"<<endl;
    of<<"NACA"<<NAt<<endl;
    of<<"�嗃�A�X�y�N�g��:"<<ARt<<endl;
    of<<"��������g�͊p"<<a0t<<"[deg]"<<endl;
    of<<"�����}�p(��)"<<alphat<<"[deg]";
    of<<endl;
    of<<"-------����-------"<<endl;
    of<<"�@�̏d�ʁF"<<m<<"[g]"<<endl;
    of<<"���q���x�F"<<v<<"[m/s]"<<endl;
    of<<"-------�݌v��-------"<<endl;
    of<<"���S���F"<<as<<endl;
    of<<"�O��g�͔z��Lt�FLw=1�F"<<x<<endl;
    of<<endl;
    of<<"-------��́��݌v����-------"<<endl;
    of<<"-------���ʐ�-------"<<endl;
    of<<"�嗃�S��"<<bw<<"[mm]"<<endl;
    of<<"�嗃�R�[�h��"<<cw<<"[mm]"<<endl;
    of<<"�����S��"<<bt<<"[mm]"<<endl;
    of<<"�����R�[�h��"<<ct<<"[mm]"<<endl;
    of<<"�嗃�ʐ�"<<Sw<<"[mm^2]"<<endl;
    of<<"�����ʐ�"<<St<<"[mm^2]"<<endl;
    of<<endl;
    of<<"-------�����g��-------"<<endl;
    of<<"�嗃���C�m���Y���F"<<Rew<<endl;
    of<<"�g�͌W���F"<<CLw<<endl;
    of<<"�嗃�����g��:"<<Lwg<<"[g]"<<endl;
    of<<"�������C�m���Y���F"<<Ret<<endl;
    of<<"�g�͌W���F"<<CLt<<endl;
    of<<"���������g��:"<<Ltg<<"[g]"<<endl;
    of<<"���v�����g��:"<<Lwg+Ltg<<"[g]"<<endl;
    of<<endl;
    of<<"-------����݌v-------"<<endl;
    of<<"CL��:"<<CLalpha<<endl;
    of<<"�݌v�d�S�ʒu"<<hn/2<<"[%]"<<endl;
    of<<"�嗃��͒��S����̏d�S�ʒu"<<lw<<"[mm]"<<endl;
    of<<"������͒��S����̏d�S�ʒu"<<lt<<"[mm]"<<endl;
    of<<"�嗃�O������̏d�S�ʒu"<<hn*cw/2<<"[mm]"<<endl;
    of<<"�����O������̏d�S�ʒu"<<lt+ct*hnw<<"[mm]"<<endl;
    of<<endl;
    of<<"-------�r���v�Z����-------"<<endl;
    of<<"�@�̏d��:"<<mkg<<"[Kg]"<<endl;
    of<<"�@�̏d�ʁF"<<mN<<"[N]"<<endl;
    of<<"�S���S�����g��"<<mas<<"[N]"<<endl;
    of<<"������͒��S����̏d�S�ʒu"<<lt<<"[mm]"<<endl;


  return(0);
}

int main(){
	  int p;
	  Airplane I;
	  I.iuputdata();
	  I.analysisdata();
	  cout<<"�v�Z���ʂ�ۑ����܂����H(Yes:1�@No:�Q)";
	  cin>>p;

	  if(p==1){
	      I.outputdata();
	  }

	  return(0);
}
