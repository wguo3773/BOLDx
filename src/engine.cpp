#include <Rcpp.h>
using namespace Rcpp;
using std::vector;

vector<double> pava(vector<double> x, vector<double> w) {
  int n=x.size(); vector<int> block(n); for(int i=0;i<n;i++)block[i]=i;
  for(int iter=0;iter<n*n;iter++) {
    int k=-1; for(int i=0;i<n-1;i++)if(x[i]>x[i+1]+1e-14){k=i;break;}
    if(k<0)break;
    int a=block[k],b=block[k+1]; double sw=0,sx=0;
    for(int i=0;i<n;i++)if(block[i]==a||block[i]==b){sw+=w[i];sx+=w[i]*x[i];}
    double z=sw>0?sx/sw:(x[k]+x[k+1])/2;
    for(int i=0;i<n;i++)if(block[i]==a||block[i]==b){x[i]=z;block[i]=a;}
  }
  return x;
}

// Sequential child-specific plug-in scores, not Bayesian marginal evidence.
// Contributions are processed in enrollment order; merge equal-count components.
vector<double> mixweights(double A,double B,vector<double> q,vector<int> lab,double se,double sp) {
  vector<double>w(1,1);
  for(unsigned i=0;i<q.size();i++) {
    double l1=lab[i]<0?1:(lab[i]?se:1-se),l0=lab[i]<0?1:(lab[i]?1-sp:sp);
    vector<double> next(w.size()+1,0);
    for(unsigned k=0;k<w.size();k++) {
      double a=A+k,b=B+i-k,den=a+b+1;
      double theta1=(a+1)/den,theta0=a/den;
      double den1=l1*q[i]*theta1+l0*(1-theta1);
      double den0=l1*q[i]*theta0+l0*(1-theta0);
      double score1=den1>0?l1*q[i]*theta1/den1:0;
      double score0=den0>0?l0*(1-theta0)/den0:0;
      next[k]+=w[k]*score0;
      next[k+1]+=w[k]*score1;
    }
    double den=0;for(double z:next)den+=z;
    if(den<=0)stop("Zero mixture normalizer");
    for(double &z:next)z/=den; w=next;
  }
  return w;
}

// [[Rcpp::export]]
List posterior_check(double A,double B,NumericVector q,IntegerVector labels,double se=.7,double sp=.9,double phi=.3) {
  auto w=mixweights(A,B,as<vector<double>>(q),as<vector<int>>(labels),se,sp);
  double mean=0,cp=0; int m=q.size();
  for(int k=0;k<=m;k++){mean+=w[k]*(A+k)/(A+B+m);cp+=w[k]*R::pbeta(phi,A+k,B+m-k,false,false);}
  return List::create(_["weights"]=w,_["mean"]=mean,_["cp"]=cp);
}

int boin_final(vector<int> n,vector<int> y,vector<bool> excluded,double phi,double gamma) {
  int J=n.size();vector<double> vals,weights;vector<int> idx;
  for(int j=0;j<J;j++)if(n[j]>=3 && R::pbeta(phi,y[j]+1,n[j]-y[j]+1,false,false)>gamma)
    for(int k=j;k<J;k++)excluded[k]=true;
  if(excluded[0])return 0;
  for(int j=0;j<J;j++)if(n[j]>0&&!excluded[j]) {
    double a=y[j]+.05,b=n[j]-y[j]+.05;
    vals.push_back(a/(a+b));weights.push_back((a+b)*(a+b)*(a+b+1)/(a*b));idx.push_back(j);
  }
  if(idx.empty())return 0;
  vals=pava(vals,weights);int best=0;double loss=1e9;
  for(unsigned k=0;k<vals.size();k++){double d=std::abs(vals[k]+(k+1)*1e-10-phi);if(d<loss){loss=d;best=k;}}
  return idx[best]+1;
}

// [[Rcpp::export]]
int final_boin_check(IntegerVector n,IntegerVector y,double phi=.3,double gamma=.95) {
  return boin_final(as<vector<int>>(n),as<vector<int>>(y),vector<bool>(n.size(),false),phi,gamma);
}

// [[Rcpp::export]]
int tite_action_check(int n,int y,int pending,double stft,double phi=.3) {
  double p1=.6*phi,p2=1.4*phi;
  double le=log((1-p1)/(1-phi))/log(phi*(1-p1)/(p1*(1-phi)));
  double ld=log((1-phi)/(1-p2))/log(p2*(1-phi)/(phi*(1-p2)));
  if(!pending)return (double)y/n<=le?1:((double)y/n>=ld?-1:0);
  double p=(y+phi/2)/(n-pending+1),est=(y+p/(1-p)*(pending-stft))/n;
  return (double)y/n<phi-1e-12&&est<=le?1:((double)y/n>=phi-1e-12&&est>=ld?-1:0);
}

// mode: 0=TITE-BOIN, 1=BOLD follow-up, 2=legacy prediction (not reported),
// 3=observed BOIN, 4=observed BOLD, 5=intervention only, 6=intervention plus follow-up only.
// All simulation truths are isolated here; posterior updates receive observable histories only.
// [[Rcpp::export]]
NumericMatrix calendar_combined(NumericVector truth,double phi,int nsim=10000,int mode=0,int seed=20261014,
 int nmax=36,double window=3,double rate=2,double gamma=.95,double tau=.49,
 double sensitivity=.7,double specificity=.9,int timing=0,int early_n=100,
 bool min_two=true,bool trace=false,double eligibility=1,double efficacy=1,
 Nullable<NumericVector> prior_means=R_NilValue,
 Nullable<NumericVector> prior_sizes=R_NilValue,int start=1,
 double pending_limit=.5,int arrival=0) {
  Function setseed("set.seed");setseed(seed);int J=truth.size();
  NumericMatrix out(nsim,19+2*J);
  CharacterVector names={"selected","n","dlt","duration","suspension","suspensions","max_pending",
    "mean_pending_decision","decisions","safety_stop","last_enrollment","above_mtd","at_mtd"};
  for(int j=0;j<J;j++)names.push_back("n"+std::to_string(j+1));
  for(int j=0;j<J;j++)names.push_back("y"+std::to_string(j+1));
  names.push_back("latent_dlt");names.push_back("prevented");names.push_back("interventions");names.push_back("eligible");names.push_back("intervention_failures");names.push_back("mixture_patients");
  colnames(out)=names;
  int true_mtd=0;double md=10;for(int j=0;j<J;j++)if(std::abs(truth[j]-phi)<md){md=std::abs(truth[j]-phi);true_mtd=j;}
  // Lookup tables depend only on prior, counts and target, never true dose risks.
  NumericVector pm=prior_means.isNotNull()?NumericVector(prior_means):NumericVector(J,phi);
  NumericVector ps=prior_sizes.isNotNull()?NumericVector(prior_sizes):NumericVector(J,3.0);
  vector<double>A(J),B(J);
  vector<vector<double>> boin_tail(nmax+1,vector<double>(nmax+1));
  vector<vector<vector<double>>> tails(J,boin_tail);
  for(int j=0;j<J;j++) {
    A[j]=pm[j]*ps[j];B[j]=(1-pm[j])*ps[j];
    for(int n=0;n<=nmax;n++)for(int y=0;y<=n;y++)
      tails[j][n][y]=R::pbeta(phi,A[j]+y,B[j]+n-y,false,false);
  }
  for(int n=0;n<=nmax;n++)for(int y=0;y<=n;y++)
    boin_tail[n][y]=R::pbeta(phi,y+1,n-y+1,false,false);
  for(int rep=0;rep<nsim;rep++) {
    if(rep%1000==0)checkUserInterrupt();
    vector<double> u(nmax),v(nmax),l(nmax),gap(nmax);
    for(int i=0;i<nmax;i++){u[i]=R::runif(0,1);v[i]=R::runif(0,1);l[i]=R::runif(0,1);gap[i]=arrival?R::rexp(1/rate):R::runif(0,2/rate);}
    vector<double>eligible_u(nmax),success_u(nmax);
    for(int i=0;i<nmax;i++){eligible_u[i]=R::runif(0,1);success_u[i]=R::runif(0,1);}
    vector<int> dose,y,label,intervened,latent,has_prediction;vector<double> enrolled,event,resolve;vector<int> n(J,0);
    int nelig=0,nint=0,nprevent=0,nfail=0;
    vector<bool> excluded(J,false);int current=start-1,N=0,suspensions=0,maxpending=0,decisions=0,pendsum=0;
    double time=0,wait=0;bool safety=false,done=false;
    while(!done) {
      // Within each cohort, patients enter sequentially at a fixed assigned dose.
      for(int c=0;c<3&&N<nmax;c++) {
        if(excluded[current])stop("Attempted enrollment at an excluded dose");
        if(c>0)time+=gap[N];
        int yi=u[N]<truth[current];double ti=window;
        if(yi) {
          if(timing==1)ti=window*v[N];
          else {
            // Weibull calibrated to P(T<=W)=p and conditional median W/2.
            double p=truth[current],a=-log(1-p),b=-log(1-p*.5);
            double shape=log(b/a)/log(.5);
            ti=window*pow(-log(1-v[N]*p)/a,1/shape);
            if(timing==2)ti=window*sqrt(v[N]); // Deliberately late conditional timing stress test.
          }
        }
        int li=l[N]<(yi?sensitivity:1-specificity);
        bool scr=(mode==5||mode==6)&&eligible_u[N]<eligibility;
        bool inter=scr&&li;int prevent=yi&&inter&&success_u[N]<efficacy;
        nelig+=scr;nint+=inter;nprevent+=prevent;
        int actual=yi-prevent; nfail+=inter&&actual;
        dose.push_back(current);y.push_back(actual);latent.push_back(yi);label.push_back(li);intervened.push_back(inter);
        has_prediction.push_back(scr);
        enrolled.push_back(time);event.push_back(actual?time+ti:INFINITY);resolve.push_back(time+(actual?ti:window));
        n[current]++;N++;
        int now_pending=0;for(int i=0;i<N;i++)if(resolve[i]>time+1e-10)now_pending++;
        maxpending=std::max(maxpending,now_pending);
      }
      if(N>=nmax)break;
      // Evaluate at the next candidate's arrival, using the latest observable data.
      time+=gap[N];
      bool chosen=false,paused=false;
      while(!chosen&&!done) {
        // Protocol assumption: intervention follow-up is resolved before another
        // dose is assigned. Wait using observable follow-up for EVERY treated
        // intervention patient, without using their latent untreated status.
        if(mode==6) {
          double intervention_resolution=time;
          for(int i=0;i<N;i++)if(intervened[i])
            intervention_resolution=std::max(intervention_resolution,resolve[i]);
          if(intervention_resolution>time+1e-10) {
            wait+=intervention_resolution-time;time=intervention_resolution;
            if(!paused)suspensions++;paused=true;
          }
        }
        vector<int>x(J,0),complete(J,0),pending(J,0);vector<double>stft(J,0),cp(J),means(J);
        vector<vector<double>>q(J);vector<vector<int>> labs(J);
        int total_pending=0;double next=INFINITY;
        for(int i=0;i<N;i++) {
          int j=dose[i];
          if(intervened[i] && event[i]>time+1e-10) {
            // Only intervention patients without an observed DLT branch.
            // Observed intervention failures enter x and complete below.
            q[j].push_back(1);labs[j].push_back(1);
            if(resolve[i]>time+1e-10){pending[j]++;total_pending++;if(j==current)next=std::min(next,resolve[i]);}
            continue;
          }
          if(resolve[i]<=time+1e-10){complete[j]++;if(event[i]<=time+1e-10)x[j]++;}
          else {pending[j]++;total_pending++;stft[j]+=std::min(1.0,(time-enrolled[i])/window);
            q[j].push_back(std::max(0.0,1-(time-enrolled[i])/window));labs[j].push_back((mode==2||(mode==6&&has_prediction[i]))?label[i]:-1);
            if(j==current)next=std::min(next,resolve[i]);}
        }
        maxpending=std::max(maxpending,total_pending);
        if((mode==3||mode==4||mode==5) && total_pending>0) {
          double last=time;for(int i=0;i<N;i++)last=std::max(last,resolve[i]);
          wait+=last-time;time=last;if(!paused)suspensions++;paused=true;continue;
        }
        bool boin=mode==0||mode==3;
        for(int j=0;j<J;j++) {
          if(boin) {cp[j]=boin_tail[n[j]][x[j]];}
          else {
            auto w=mixweights(A[j]+x[j],B[j]+complete[j]-x[j],q[j],labs[j],sensitivity,specificity);
            cp[j]=0;means[j]=0;
            for(unsigned k=0;k<w.size();k++){cp[j]+=w[k]*tails[j][n[j]][x[j]+k];means[j]+=w[k]*(A[j]+x[j]+k)/(A[j]+B[j]+n[j]);}
          }
          if(n[j]>=3 && cp[j]>(boin?gamma:gamma-1e-12))for(int k=j;k<J;k++)excluded[k]=true;
        }
        if(excluded[0]){safety=true;done=true;break;}
        int candidate=current,action=0;
        if(boin) {
          action=tite_action_check(n[current],x[current],pending[current],stft[current],phi);
          if(excluded[current])action=-1;
          if(action==1&&current<J-1&&!excluded[current+1])candidate=current+1;
          if(action==-1&&current>0)candidate=current-1;
        } else {
          vector<double>local=cp,weights(n.begin(),n.end());
          for(int j=0;j<J;j++){if(j<current-1)local[j]=0;if(j>current+1)local[j]=100;}
          auto pp=pava(local,weights);for(double &z:pp)z=std::round(z*1000)/1000;
          double loss=1e9;vector<int> ties;
          for(int j=std::max(0,current-1);j<=std::min(J-1,current+1);j++)if(!excluded[j]) {
            double d=std::abs(pp[j]-tau);if(d<loss-1e-12){loss=d;ties={j};}else if(std::abs(d-loss)<1e-12)ties.push_back(j);
          }
          if(!ties.empty())candidate=pp[ties[0]]<tau?ties.back():ties.front();
        }
        // Late events at an earlier dose can eliminate several levels at once.
        if(excluded[candidate]) {
          candidate=0;for(int j=0;j<J;j++)if(!excluded[j])candidate=j;
        }
        bool observed=mode==3||mode==4||mode==5;
        bool suspend=observed?total_pending>0:pending[current]>pending_limit*n[current];
        // The optional guard counts actual completed follow-up, including intervention.
        if(min_two&&candidate>current&&n[current]-pending[current]<2)suspend=true;
        // Safety de-escalation/elimination does not require waiting for outcomes.
        double ld=log((1-phi)/(1-1.4*phi))/log(1.4*(1-phi)/(1-1.4*phi));
        if(!observed&&(excluded[current]||(boin&&x[current]/(double)n[current]>=ld)))suspend=false;
        if(suspend) {
          if(observed){next=INFINITY;for(int i=0;i<N;i++)if(resolve[i]>time+1e-10)next=std::min(next,resolve[i]);}
          if(!R_finite(next))stop("No future event to release suspension");
          wait+=next-time;time=next;if(!paused)suspensions++;paused=true;continue;
        }
        decisions++;pendsum+=total_pending;
        if(trace)Rcout<<"rep "<<rep<<" time "<<time<<" N "<<N<<" current "<<current+1<<" pending "<<total_pending<<" next "<<candidate+1<<"\n";
        if(n[current]>=early_n&&candidate==current&&!excluded[current]){done=true;break;}
        current=candidate;chosen=true;
      }
    }
    vector<int>x(J,0);for(int i=0;i<N;i++)x[dose[i]]+=y[i];
    double duration=time;for(int i=0;i<N;i++)duration=std::max(duration,resolve[i]);
    int selected=0;
    if(!safety) {
      if(mode==0||mode==3)selected=boin_final(n,x,excluded,phi,gamma);
      else {
        vector<double>means(J),weights(n.begin(),n.end());
        for(int j=0;j<J;j++) {
          int o=0,xx=0,m=0;
          for(int i=0;i<N;i++)if(dose[i]==j){if(intervened[i] && !y[i])m++;else{o++;xx+=y[i];}}
          auto w=mixweights(A[j]+xx,B[j]+o-xx,vector<double>(m,1),vector<int>(m,1),sensitivity,specificity);
          double cp=0;means[j]=0;
          for(int k=0;k<=m;k++){means[j]+=w[k]*(A[j]+xx+k)/(A[j]+B[j]+n[j]);cp+=w[k]*tails[j][n[j]][xx+k];}
          if(n[j]>=3&&cp>=gamma)for(int k=j;k<J;k++)excluded[k]=true;
        }
        if(!excluded[0]) {
          auto fit=pava(means,weights);for(double&z:fit)z=std::round(z*1000)/1000;
          double loss=1e9;vector<int>ties;
          for(int j=std::max(0,current-1);j<=std::min(J-1,current+1);j++)if(n[j]>0&&!excluded[j]) {
            double d=std::abs(fit[j]-phi);if(d<loss-1e-12){loss=d;ties={j};}else if(std::abs(d-loss)<1e-12)ties.push_back(j);
          }
          if(!ties.empty())selected=(fit[ties.front()]<phi&&fit[ties.back()]<phi?ties.back():ties.front())+1;
        }
      }
    }
    int ndlt=0,over=0;for(int j=0;j<J;j++){ndlt+=x[j];if(j>true_mtd)over+=n[j];}
    out(rep,0)=selected;out(rep,1)=N;out(rep,2)=ndlt;out(rep,3)=duration;
    out(rep,4)=wait;out(rep,5)=suspensions;out(rep,6)=maxpending;
    out(rep,7)=decisions?(double)pendsum/decisions:0;out(rep,8)=decisions;out(rep,9)=safety;
    out(rep,10)=enrolled.back();out(rep,11)=over;out(rep,12)=n[true_mtd];
    for(int j=0;j<J;j++){out(rep,13+j)=n[j];out(rep,13+J+j)=x[j];}
    int latent_total=0;for(int z:latent)latent_total+=z;
    if(ndlt+nprevent!=latent_total)stop("DLT accounting failure");
    if(nfail>ndlt || nfail+nprevent>nint)stop("Intervention accounting failure");
    out(rep,13+2*J)=latent_total;out(rep,14+2*J)=nprevent;out(rep,15+2*J)=nint;out(rep,16+2*J)=nelig;out(rep,17+2*J)=nfail;out(rep,18+2*J)=nint-nfail;
  }
  return out;
}
