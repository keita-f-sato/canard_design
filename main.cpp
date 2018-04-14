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
    cout<<"-------翼最適化設計-------"<<endl;
    cout<<"-------主翼-------"<<endl;
    cout<<"主翼の各データを入力してください"<<endl;
    cout<<"NACA";
    cin>>NAw;
    cout<<"主翼アスペクト比:";
    cin>>ARw;
    cout<<"主翼無揚力角[deg]:";
    cin>>a0w;
    cout<<"主翼迎角(仮)[deg]";
    cin>>alphaw;
    cout<<endl;

    cout<<"-------先尾翼-------"<<endl;
    cout<<"先尾翼の各データを入力してください"<<endl;
    cout<<"NACA";
    cin>>NAt;

    cout<<"先尾翼アスペクト比:";
    cin>>ARt;
    cout<<"先尾翼無揚力角[deg]:";
    cin>>a0t;
    cout<<"尾翼迎角(仮)[deg]";
    cin>>alphat;
    cout<<endl;

    cout<<"-------諸元-------"<<endl;
    cout<<"機体諸元を入力してください"<<endl;
    cout<<"機体重量[g]：";
    cin>>m;
    cout<<"巡航速度[m/s]：";
    cin>>v;
    cout<<"-------設計率-------"<<endl;
    cout<<"安全率：";
    cin>>as;
    cout<<"前後揚力配分Lt：Lw=1：";
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

    cout<<"-------解析＆設計結果-------"<<endl;
    cout<<"-------翼面積-------"<<endl;
    cout<<"主翼全長"<<bw<<endl;
    cout<<"主翼コード長"<<cw<<endl;
    cout<<"主翼面積"<<Sw<<endl;
    cout<<"尾翼全長"<<bt<<endl;
    cout<<"尾翼コード長"<<ct<<endl;
    cout<<"尾翼面積"<<St<<endl;
    cout<<"-------発生揚力-------"<<endl;
    cout<<"主翼レイノルズ数[-]"<<Rew<<endl;
    cout<<"主翼揚力係数[-]"<<CLw<<endl;
    cout<<"主翼発生揚力[g]"<<Lwg<<endl;
    cout<<"主翼揚力係数[-]"<<CLw<<endl;
    cout<<"尾翼レイノルズ数[-]"<<Ret<<endl;
    cout<<"尾翼揚力係数[-]"<<CLt<<endl;
    cout<<"尾翼発生揚力[g]"<<Ltg<<endl;
    cout<<"合計発生揚力[g]"<<Lwg+Ltg<<endl;
    cout<<"-------安定設計-------"<<endl;
    cout<<"CLα:"<<CLalpha<<endl;
    cout<<"設計重心位置[%]"<<hn/2<<endl;
    cout<<"主翼空力中心からの重心位置[mm]"<<lw<<endl;
    cout<<"尾翼空力中心からの重心位置[mm]"<<lt1<<endl;
    cout<<"主翼前縁からの重心位置[mm]"<<cw*(hn/2)<<endl;
    cout<<"尾翼前縁からの重心位置[mm]"<<lt1+ct*hnw<<endl;

   return(0);
}

int Airplane::outputdata(){

    ofstream of;
    of.open("Airplane.txt",ios::app);
    of<<endl;
    of<<"-------翼最適化設計-------"<<endl;
    of<<"-------主翼-------"<<endl;
    of<<"NACA"<<NAw<<endl;
    of<<"主翼アスペクト比:"<<ARw<<endl;
    of<<"主翼無揚力角"<<a0w<<"[deg]"<<endl;
    of<<"主翼迎角(仮)"<<alphaw<<"[deg]"<<endl;
    of<<endl;
    of<<"-------先尾翼-------"<<endl;
    of<<"NACA"<<NAt<<endl;
    of<<"主翼アスペクト比:"<<ARt<<endl;
    of<<"先尾翼無揚力角"<<a0t<<"[deg]"<<endl;
    of<<"尾翼迎角(仮)"<<alphat<<"[deg]";
    of<<endl;
    of<<"-------諸元-------"<<endl;
    of<<"機体重量："<<m<<"[g]"<<endl;
    of<<"巡航速度："<<v<<"[m/s]"<<endl;
    of<<"-------設計率-------"<<endl;
    of<<"安全率："<<as<<endl;
    of<<"前後揚力配分Lt：Lw=1："<<x<<endl;
    of<<endl;
    of<<"-------解析＆設計結果-------"<<endl;
    of<<"-------翼面積-------"<<endl;
    of<<"主翼全長"<<bw<<"[mm]"<<endl;
    of<<"主翼コード長"<<cw<<"[mm]"<<endl;
    of<<"尾翼全長"<<bt<<"[mm]"<<endl;
    of<<"尾翼コード長"<<ct<<"[mm]"<<endl;
    of<<"主翼面積"<<Sw<<"[mm^2]"<<endl;
    of<<"尾翼面積"<<St<<"[mm^2]"<<endl;
    of<<endl;
    of<<"-------発生揚力-------"<<endl;
    of<<"主翼レイノルズ数："<<Rew<<endl;
    of<<"揚力係数："<<CLw<<endl;
    of<<"主翼発生揚力:"<<Lwg<<"[g]"<<endl;
    of<<"尾翼レイノルズ数："<<Ret<<endl;
    of<<"揚力係数："<<CLt<<endl;
    of<<"尾翼発生揚力:"<<Ltg<<"[g]"<<endl;
    of<<"合計発生揚力:"<<Lwg+Ltg<<"[g]"<<endl;
    of<<endl;
    of<<"-------安定設計-------"<<endl;
    of<<"CLα:"<<CLalpha<<endl;
    of<<"設計重心位置"<<hn/2<<"[%]"<<endl;
    of<<"主翼空力中心からの重心位置"<<lw<<"[mm]"<<endl;
    of<<"尾翼空力中心からの重心位置"<<lt<<"[mm]"<<endl;
    of<<"主翼前縁からの重心位置"<<hn*cw/2<<"[mm]"<<endl;
    of<<"尾翼前縁からの重心位置"<<lt+ct*hnw<<"[mm]"<<endl;
    of<<endl;
    of<<"-------途中計算結果-------"<<endl;
    of<<"機体重量:"<<mkg<<"[Kg]"<<endl;
    of<<"機体重量："<<mN<<"[N]"<<endl;
    of<<"全安全発生揚力"<<mas<<"[N]"<<endl;
    of<<"尾翼空力中心からの重心位置"<<lt<<"[mm]"<<endl;


  return(0);
}

int main(){
	  int p;
	  Airplane I;
	  I.iuputdata();
	  I.analysisdata();
	  cout<<"計算結果を保存しますか？(Yes:1　No:２)";
	  cin>>p;

	  if(p==1){
	      I.outputdata();
	  }

	  return(0);
}
