### Return Bias-Corrected Log Periodograms
LogPdg = function(inputData){
  data1 = as.matrix(inputData)
  
  ch = dim(data1)[1];
  T = dim(data1)[2];
  J = floor((T+1)/2)-1;
  
  ##Calculate log periodogram
  logp = rep(0,ch*J);
  for(i in 1:ch) {
    logp[((i-1)*J+1):(i*J)] = log((1/T)*(abs(fft(data1[i,])[1:J]))^2);
  }
  
  return(logp+0.57721);
}

### Periodogram not log periodograms
RawPdg = function(inputData){
  data1 = as.matrix(inputData)
  ch=dim(data1)[1];
  T=dim(data1)[2];
  J = floor((T+1)/2)-1;
  Rawp = rep(0,ch*J);
  for(i in 1:ch) {
    Rawp[((i-1)*J+1):(i*J)] =(1/T)*(abs(fft(data1[i,])[1:J]))^2
  }
  return(Rawp);
}

### Convert Vector (Log Periodogram) to Matrix
### Columns are objects, with time series are stacked in rows.
VecToMatrix = function(inputData, ch, T, J){
  matrixp=matrix(0,T/2-1,ch)
  J=T/2-1
  for(i in 1:ch){
    matrixp[,i]= inputData[((i-1)*J+1):(i*J)]
  }
  return(matrixp)
}

#define neighborhood
neighborhood = function(x,x.star,y,bandwidth){
  index=c(1:length(x))
  indexn = index[(x<(x.star+bandwidth+1)) & (x>(x.star-bandwidth-1))]
  xn = x[indexn]
  yn = y[indexn]
  out=list(xn,yn)
  names(out) = c("xn", "yn")
  return(out)
}

#Calculate local estimate
localEstimates=function(xn,yn,bandwith){
  #boxcar smoother W = 1/(2p+1)
  w=rep(1/(2*bandwith + 1),length(xn))
  f.hat = as.numeric((t(w)%*%yn)/sum(w) )
  #f.hat = sum(yn)/(2*bandwith + 1)
  return(f.hat)
}

### Average Channel by using local channel around it.
AvgChannel = function(inputData,channel){
  for(i in 1:160){
    inputData[i,channel] = mean(c(inputData[i,(channel-5):(channel-1)],inputData[i,(channel+1):(channel+5)]))
  }
  return(inputData)
}

SmoothGVC = function(matrix.Rawpdg,span){
  GVCp = rep(0,length(span))
  span.min = 0
  n.trial = dim(matrix.Rawpdg)[1]
  matrix.smooth = matrix(0,dim(matrix.Rawpdg)[1],dim(matrix.Rawpdg)[2])
  y.estimate = matrix(0,length(span),dim(matrix.Rawpdg)[2])
  for(k in 1:n.trial){
    for(i in 1:(length(span))){
      bandwidth = span[i]
      first.trial.process = ProcessRawPdg(matrix.Rawpdg[k,],bandwidth)
      x.feq = seq(1,length(first.trial.process),1)
      for(j in 1:length(matrix.Rawpdg[k,])){
        window = neighborhood(x=x.feq,x.star=x.feq[j]+bandwidth,y=first.trial.process,bandwidth)
        y.estimate[i,j]=localEstimates(window$xn,window$yn,bandwidth)
      }
      GVCp[i] = CalculateGVC(matrix.Rawpdg[k,], y.estimate[i,], bandwidth)
    }
    matrix.smooth[k,] = y.estimate[which(GVCp == min(GVCp)),]
  }  
  return(matrix.smooth)
}

CalculateGVC = function(f, f.hat, bandwidth){
  M = length(f)
  sum = 0
  q = c(0.5,rep(1,M-2),0.5)
  for(i in 1:M){
    num = -log(f[i]/f.hat[i])+(f[i]-f.hat[i])/f.hat[i]
    dem = (1 - (1/(2*bandwidth + 1)))^2
    sum = sum + q[i]*(num/dem)
  }
  return(sum/M)
}

ProcessRawPdg = function(pdg,bandwidth){
  temp = rep(0,bandwidth)
  end = length(pdg)
  temp = rev(pdg[2:(bandwidth+1)])
  pdg.final = c(temp,pdg)
  temp = rev(pdg[(end-bandwidth-1):(end-1)])
  pdg.final = c(pdg.final,temp)
  return(pdg.final)
}

LogSmoothGVC = function(inputdata,span,T,J,ch){
  Channel.RawPdg = RawPdg(inputdata)
  Channel.RawPdg=t(VecToMatrix(Channel.RawPdg,ch,T,J))
  #Channel.RawPdg=AvgChannel(Channel.RawPdg,61,ch) #filter out 61 Hz, not needed
  return(log(SmoothGVC(Channel.RawPdg,span)) + 0.57721)
}

######################contamination models########
generatesim=function(type=1, rate=0.2, T=400, neach=150, amp=8){
  #type1:eyeblink
  #type2:eyemovement
  ntrial = neach/5
  ch = neach*4
  nch = 20
  J = floor((T+1)/2)-1
  a = rnorm(5,0,0.001)#random difference
  a1 = a[1]
  a2 = a[2]
  a3 = a[3]
  a4 = a[4]
  a5 = a[5]
  data1 = matrix(nrow=neach, ncol=T);   #cluster 1 
  for(i in 1:(1*neach/5)){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a1,0.1)),n=T)}
  for(i in (1*neach/5+1):(2*neach/5)){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a2,0.1)),n=T)}
  for(i in (2*neach/5+1):(3*neach/5)){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a3,0.1)),n=T)}
  for(i in (3*neach/5+1):(4*neach/5)){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a4,0.1)),n=T)}
  for(i in (4*neach/5+1):(5*neach/5)){
    data1[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.8+a5,0.1)),n=T)}
  
  data2 = matrix(nrow=neach,ncol=T)
  for(i in 1:(1*neach/5)){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a1,-0.9)),n=T)}
  for(i in (1*neach/5+1):(2*neach/5)){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a2,-0.9)),n=T)}
  for(i in (2*neach/5+1):(3*neach/5)){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a3,-0.9)),n=T)}
  for(i in (3*neach/5+1):(4*neach/5)){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a4,-0.9)),n=T)}
  for(i in (4*neach/5+1):(5*neach/5)){
    data2[i,] = arima.sim(list(order=c(2,0,0),ar=c(0.9+a5,-0.9)),n=T)}
  
  
  data3 = matrix(nrow=neach,ncol=T)
  for(i in 1:(1*neach/5)){
    data3[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.1+a1,-0.9)),n=T)}
  for(i in (1*neach/5+1):(2*neach/5)){
    data3[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.1+a2,-0.9)),n=T)}
  for(i in (2*neach/5+1):(3*neach/5)){
    data3[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.1+a3,-0.9)),n=T)}
  for(i in (3*neach/5+1):(4*neach/5)){
    data3[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.1+a4,-0.9)),n=T)}
  for(i in (4*neach/5+1):(5*neach/5)){
    data3[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.1+a5,-0.9)),n=T)}
  
  
  data4 = matrix(nrow=neach,ncol=T)
  for(i in 1:(1*neach/5)){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a1,-0.9)),n=T)}
  for(i in (1*neach/5+1):(2*neach/5)){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a2,-0.9)),n=T)}
  for(i in (2*neach/5+1):(3*neach/5)){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a3,-0.9)),n=T)}
  for(i in (3*neach/5+1):(4*neach/5)){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a4,-0.9)),n=T)}
  for(i in (4*neach/5+1):(5*neach/5)){
    data4[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.9+a5,-0.9)),n=T)}
  
  
  data5 = matrix(nrow=neach,ncol=T)
  for(i in 1:(1*neach/5)){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a1,0.1)),n=T)}
  for(i in (1*neach/5+1):(2*neach/5)){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a2,0.1)),n=T)}
  for(i in (2*neach/5+1):(3*neach/5)){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a3,0.1)),n=T)}
  for(i in (3*neach/5+1):(4*neach/5)){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a4,0.1)),n=T)}
  for(i in (4*neach/5+1):(5*neach/5)){
    data5[i,] = arima.sim(list(order=c(2,0,0),ar=c(-0.8+a5,0.1)),n=T)}
  
  
  data11 = data1
  data22 = (4/5)*data1+(1/10)*data2
  data33 = (3/5)*data1+(1/10)*data3
  data44 = (2/5)*data1+(1/10)*data4
  data55 = (1/5)*data1+(1/10)*data5
  data=rbind(data11,data22,data33,data44)
  
  span = seq(from = 3, to = 15, by = 1)
  if (type==0){
    simm = LogSmoothGVC(inputdata =data,span,T,J,ch)
    sim = array(0,c(ntrial,J,nch))
    for(i in 1:nch){
      sim[,,i] = simm[((i-1)*ntrial+1):(i*ntrial),]
    }
  }
  
  if(type==1){

    dataout =  eyeblk(data, rate, amp)
    simmout = LogSmoothGVC(inputdata =dataout,span,T,J,ch)
    sim = array(0,c(ntrial,J,nch))
    for(i in 1:nch){
      sim[,,i] = simmout[((i-1)*ntrial+1):(i*ntrial),]
    }
  }
  if(type==2){
    dataout =  eyemvt(data, rate, 16)
    simmout = LogSmoothGVC(inputdata =dataout,span,T,J,ch)
    sim = array(0,c(ntrial,J,nch))
    for(i in 1:nch){
      sim[,,i] = simmout[((i-1)*ntrial+1):(i*ntrial),]
    }
  }
  
  return(sim)
}

#generate data in Experiment1
generatesim2=function(type=1, rate=0.2, T=200, neach=150){
  nclu = 4
  ntrial = neach/5
  ch = neach*4
  nch = 20
  t = (1:T)/T
  a = rnorm(5,0,0.01)#random difference
  
  
  sim = array(0,c(ntrial,T,nch))
  for(i in 1:nclu){
    for(j in 1:(nch/nclu)){
      pro2 = runif(1)    # up shift or down shift
      for(k in 1:ntrial){
        pro = runif(1) #comtamination rate
        if(type==0){
          noise = gaussprocess(T, K=function(s, t) { exp(-abs(t-s)) })
        }
        
        if(type==1){
          if(pro<=rate){
            noise = gaussprocess(T, K=function(s, t) { exp(-abs(t-s)) }) + 8*sgn(pro2) 
          }
          if(pro>rate){
            noise = gaussprocess(T, K=function(s, t) { exp(-abs(t-s)) }) 
          }
        }
        if(type==2){
          if(pro<=rate){
            #noise1 = gaussprocess(T, K=function(s, t) { 5 *exp(-2*abs(t-s)^0.5) })
            noise1 = gaussprocess(T, K=function(s, t) { exp(-abs(t-s)) })
            noise2 = 60*(1-t)^1.5*t * sgn(pro2) 
            noise = noise2 + noise1
          }
          if(pro>rate){
            noise = gaussprocess(T, K=function(s, t) { exp(-abs(t-s)) })
          }
        }
        if(type==3){
          if(pro<=rate){
            #noise1 = gaussprocess(T, K=function(s, t) { 5 *exp(-2*abs(t-s)^0.5) })
            noise = gaussprocess(T, K=function(s, t) { 12*exp(-2*abs(t-s) ^0.1) }) + 0*sin(6*pi*t)* sgn(pro2)
          }
          if(pro>rate){
            noise = gaussprocess(T, K=function(s, t) { exp(-abs(t-s)) })
          }
        }
        sim[k, ,((i-1)*(nch/nclu)+j)] = i+2*i*t + noise
        #sim[k, ,((i-1)*(nch/nclu)+j)] = i+6*i*t + noise
      }
    }
  }
  if(type==3){
    span = seq(from = 3, to = 15, by = 1)
    J = floor((T+1)/2)-1
    
    temp = sim[,1:J,]
    for(i in 1:nch){
      temp[,,i] = LogSmoothGVC(inputdata = sim[,,i],span,T,J,ntrial)
    }
    sim = temp
  }
  return(sim)
}

#add shift to the data
addcon=function(x,rate){
  dimension=dim(x)
  m = dimension[1]
  n = dimension[2] 
  ran_v = runif(m)
  noise = rep(0,m)
  contain_num = floor(m*rate)
  index = sample(1:m, contain_num)
  for(i in  1:length(index)){
    x[index[i],]=x[index[i],]+4
  }
  return(x)
}

#add eyeblink to the data
eyeblk = function(x, rate, amp){
  dimension=dim(x)
  nch = dimension[1]
  len = dimension[2]
  
  c = floor(0.031*len):(floor(0.27*len))/(0.15*len)
  ga = gamma(c)-gamma(1.8)
  ga_1 = 2*amp*ga[length(ga):1]
  ga_2 = -amp*ga
  ga_3 = c(ga_1,ga_2)
  inter = seq(from=ga_3[floor(length(ga_3)/2)], to = ga_3[floor(length(ga_3)/2)+1], 
              by =  (ga_3[floor(length(ga_3)/2)+1] - ga_3[floor(length(ga_3)/2)]) / floor(0.032*len))
  
  gam = c()
  gam[1:length(ga_1)] = ga_1
  gam[(length(ga_1)+1):(length(ga_1)+length(inter)-2)] = inter[2:(length(inter)-1)]
  gam[(length(ga_1)+length(inter)-1):((length(ga_1)+length(inter))+length(ga_2)-2)] = ga_2
  
  contain_num = floor(nch*rate)
  index = sample(1:nch, contain_num)
  for(i in 1:length(index)){
    wn = arima.sim(list(order=c(0,0,0)),n=length(gam))
    eye = gam + wn
    x[index[i],(floor(0.25*len)+1):(floor(0.25*len)+length(gam))] = 
      x[index[i],(floor(0.25*len)+1):(floor(0.25*len)+length(gam))] + eye
  }
  
  return(x)
}

#add eyeblink to the data
eyemvt = function(x, rate, amp){
  dimension = dim(x)
  nch = dimension[1]
  len = dimension[2]
  
  index = sample(1:nch, floor(rate * nch))
  for(i in index){
    up = 1:floor(0.03 * len) * amp
    mid = rep(up[length(up)], floor(runif(1,0.1,1)*len*0.4))
    down = up[length(up):1]
    toge = c(up, mid, down)
    
    k = runif(1, 0.1,0.4)
    ind = floor(k*0.5*len):(floor(k*0.5*len) + length(toge)-1)
    x[i, ind] = x[i, ind] + toge
  }
  return(x)
}
########################clustering##################
##calculate the  central region of BD
central = function(x, tau){
  dimension = dim(x)
  m = dimension[1]
  n = dimension[2]
  depth = fMBD(x)
  
  depth_tau = quantile(depth, 1-tau)
  index = which(depth>depth_tau)
  num = x[,index]
  
  area = rep(0,m)
  for (i in 1:m){
    area[i] = max(num[i,])-min(num[i,])
  }
  return(sum(area))
}

##calculate the  central region of MS
central2 = function(x, tau){
  res_msplt = msplot1(x, plot=F)
  mo = res_msplt$mean_outlyingness
  vo = res_msplt$var_outlyingness
  out = res_msplt$outliers
  ns = dim(x)[1]-length(out)

  gs = 10
  
  if(length(out)==0 ){
    success = try('2'+'1', silent = T)
    while(class(success) == 'try-error'){
      success = try({res_c2d = c2d(mo,vo, grid_size=gs,tau, plot=F)})
      gs = gs - 1
    }
    central_curves = x[res_c2d$pts_in_pol,]
  }
  
  if(length(out)>0){
    out2 = c(out, out+length(mo)) 
    mo2 = c(mo, mo)
    vo2 = c(vo,-vo)
    x2 = rbind(x, x)
    
    success = try('2'+'1', silent = T)
    while(class(success) == 'try-error'){
      success = try({res_c2d = c2d(mo2[-out2],vo2[-out2], grid_size=gs, tau, plot=F)})
      gs = gs -1
    }
    
    central_curves = x2[-out2,][res_c2d$pts_in_pol,]
    central_curves = central_curves[1:floor(dim(central_curves)[1]/2),]
  }
  return(bandarea(central_curves))
}

##calculate the  central region of MS (SIMPLE)
central3 = function(x, tau){
  res_msplt = msplot1(x, plot=F)
  mo = res_msplt$mean_outlyingness
  vo = res_msplt$var_outlyingness
  lst = c()
  for(i in 1:length(mo)){
    if(mo[i] >= -sd(mo) & mo[i] <= sd(mo) & vo[i] <= min(vo)+sd(vo)){
      lst = append(lst,i)
    }
  }
  central_curves = x[lst,]
  return(bandarea(central_curves))
}

#area of the band
bandarea = function(a){
  len = dim(a)[2]
  are = rep(0,len)
  for(i in 1:len){
    are[i] = max(a[,i])-min(a[,i])
  }
  return(sum(are))
}

##BD distance of 2 clusters
ddcr = function(a, b,  tau, pw1=1, pw2=1){
  merg = cbind(t(a) ,t(b))
  num1 = dim(a)[1]
  num2 = dim(b)[1]
  return( (num1+num2)^pw2*central(merg, tau)^pw1 -  (num1)^pw2*central(t(a), tau)^pw1 -  (num2)^pw2*central(t(b), tau)^pw1)
}

##MS distance of 2 clusters
ddcr2 = function(a, b, tau, pw1=1, pw2=1){
  num1 = dim(a)[1]
  num2 = dim(b)[1]
  len = dim(a)[2]
  if(num1>num2){
    sam = b[sample(1:num2, num1-num2,replace = T),] + matrix(rnorm((num1-num2)*len,mean=0,sd=0.001), num1-num2,len)
    merg = cbind(t(a),t(b),t(sam))
  }
  
  if(num1<num2){
    sam = a[sample(1:num1, num2-num1,replace = T),]+ matrix(rnorm((num2-num1)*len,mean=0,sd=0.001), num2-num1,len)
    merg = cbind(t(a),t(b),t(sam))
  }
  
  
  if(num1==num2){
  merg = cbind(t(a) ,t(b))}
  merg = t(merg)
  
  success = try('2'+'1',silent = T)
  while(class(success) == 'try-error'){
    success = try({res = (num1+num2)^pw2*central2(merg, tau)^pw1 - num1^pw2*central2(a, tau)^pw1 - num2^pw1*central2(b, tau)^pw1})
    tau = tau-0.03}
  
  return(res)
}

##MS distance of 2 clusters (SIMPLE)
ddcr3 = function(a, b, tau, pw1=1, pw2=1){
  merg = cbind(t(a) ,t(b))
  merg = t(merg)
  return( (num1+num2)^pw2*central3(merg, tau)^pw1 -  (num1)^pw2*central3(t(a), tau)^pw1 -  (num2)^pw2*central3(t(b), tau)^pw1)
  return(res)
}

# DISTANCE MATRIX OF BD
distcr = function(x,tau){
  dimension = dim(x)
  p = dimension[1]
  q = dimension[2]
  r = dimension[3]
  dcr = matrix(0,r,r) 
  if(p > 1){
    for (i in 1:(r-1)){
      for(j in (i+1):r){
        dcr[i,j] = ddcr(x[,,i],x[,,j],tau)
        dcr[j,i] = dcr[i,j]####
      }
    }}
  if(p==1){
    dcr = as.matrix(dist(t(x[1,,])))
  }
  for(i in 1:r){dcr[i,i] = 10000000}  
  return(dcr)
}

# DISTANCE MATRIX OF MS
distcr2 = function(x,tau){
  dimension = dim(x)
  p = dimension[1]
  q = dimension[2]
  r = dimension[3]
  dcr = matrix(0,r,r)
  if(p > 1){
    for (i in 1:(r-1)){
      for(j in (i+1):r){
        dcr[i,j] = ddcr2(x[,,i],x[,,j],tau)
        dcr[j,i] = dcr[i,j]####
      }
    }}
  if(p==1){
    dcr = as.matrix(dist(t(x[1,,])))
  }
  for(i in 1:r){dcr[i,i] = 10000000}  
  return(dcr)
}

# DISTANCE MATRIX OF MS (SIMPLE)
distcr3 = function(x,tau){
  dimension = dim(x)
  p = dimension[1]
  q = dimension[2]
  r = dimension[3]
  dcr = matrix(0,r,r)
  if(p > 1){
    for (i in 1:(r-1)){
      for(j in (i+1):r){
        dcr[i,j] = ddcr3(x[,,i],x[,,j],tau)
        dcr[j,i] = dcr[i,j]####
      }
    }}
  if(p==1){
    dcr = as.matrix(dist(t(x[1,,])))
  }
  for(i in 1:r){dcr[i,i] = 10000000}  
  return(dcr)
}

##how many columns that is not empty
notzero = function(x){
  dimension = dim(x)
  m = dimension[1]
  n = dimension[2]
  for(i in 1:m){
    if(sum(x[i,]==rep(0,n))>5){
      z = i-1 ;break
    }
  }
  if(i==m){z=i}
  return(z)
}

###get the functional median curve
fmed = function(x){
  dep = fMBD(t(x))
  index = which.max(dep)
  s = x[index,]
  return(s)
}

#CLUSTERING
robustcluster = function(x, n, method,tau){
  dimension = dim(x)
  n_rep = dimension[1]
  n_freq = dimension[2]
  n_clu = dimension[3]
  min_dis = rep(0, n_clu)
  res = 1:n_clu
  n_element_in_each_clu = matrix(1, 1, n_clu)
  clu_mtx = array(0,c(n_clu*n_rep, n_freq,n_clu))
  for (i in 1:n_clu){
    clu_mtx[1:n_rep, , i] = x[, , i]
  }
  if(method == 'cr'){dis = distcr(x,tau)}
  if(method == 'ms'){dis = distcr2(x,tau)}
  if(method == 'ms2'){dis = distcr3(x,tau)}
  if(method == 'tv'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:(n_clu-1)){
      for(j in (i+1):n_clu){
        dis[i,j]=tvd(colMeans(x[,,i]),colMeans(x[,,j]))
        dis[j,i]=dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  if(method == 'cid'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:(n_clu-1)){
      for(j in (i+1):n_clu){
        dis[i,j]=diss.CID(colMeans(x[,,i]),colMeans(x[,,j]))
        dis[j,i]=dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  
  if(method == 'fm'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:(n_clu-1)){
      for(j in (i+1):n_clu){
        dis[i,j]=diss(colMeans(x[,,i]),colMeans(x[,,j]))
        dis[j,i]=dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  if(method == 'mean'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:n_clu){
      for(j in 1:n_clu){
        ms = rbind(colMeans(x[,,i]), colMeans(x[,,j]))
        dis[i,j] = as.matrix(dist(ms))[1,2]
        dis[j,i]=dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  if(method == 'ward'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:n_clu){
      for(j in 1:n_clu){
        ms = rbind(colMeans(x[,,i]), colMeans(x[,,j]))
        dis[i,j] = as.matrix(dist(ms))[1,2] * sqrt(n_rep^2/(2*n_rep))
        dis[j,i] = dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  
  
  rcd = 1
  for(k in n_clu:(n+1)){
    index_min = which(dis==min(dis), arr.ind=TRUE)
    index_min = index_min[1,]
    if(index_min[2] <= index_min[1]){
      temp = index_min[1]
      index_min[1] = index_min[2]
      index_min[2] = temp}
    min_dis[rcd] = dis[index_min[1], index_min[2]]
    rcd = rcd + 1
    clu_mtx[(n_rep*n_element_in_each_clu[index_min[1]]+1):(n_rep*(n_element_in_each_clu[index_min[1]]+n_element_in_each_clu[index_min[2]])), , index_min[1]]=clu_mtx[1:(n_rep*n_element_in_each_clu[index_min[2]]), , index_min[2]]
    clu_mtx[ , ,index_min[2]] = matrix(0,n_clu*n_rep, n_freq)
    n_element_in_each_clu[index_min[1]] = n_element_in_each_clu[index_min[1]]+n_element_in_each_clu[index_min[2]]
    n_element_in_each_clu[index_min[2]] = 0
    index = which(res==res[index_min[2]])
    res[index] = res[index_min[1]]
    for(i in 1:(n_clu-1)){
      for (j in (i+1):n_clu){
        a = notzero(clu_mtx[, , i])
        b = notzero(clu_mtx[, , j])
        if (a==0 || b==0){
          dis[i,j] = 10000000
          dis[j,i] = 10000000}
        else{
          if(i==index_min[1]){
            if(method == 'cr'){
              dis[i,j] = ddcr(clu_mtx[1:a, , i], clu_mtx[1:b, , j], tau)
              dis[j,i] = dis[i,j]}
            if(method == 'ms'){
              dis[i,j] = ddcr2(clu_mtx[1:a, , i], clu_mtx[1:b, , j], tau)
              dis[j,i] = dis[i,j]}
            if(method == 'ms2'){
              dis[i,j] = ddcr3(clu_mtx[1:a, , i], clu_mtx[1:b, , j], tau)
              dis[j,i] = dis[i,j]}
            if(method == 'fm'){
              dis[i,j]=diss(fmed(clu_mtx[1:a,,i]), fmed(clu_mtx[1:b,,j]))
              dis[j,i] = dis[i,j]}
            if(method == 'mean'){
              mss = rbind(colMeans(clu_mtx[1:a,,i]), colMeans(clu_mtx[1:b,,j]))
              dis[i,j] = as.matrix(dist(mss))[1,2]
              dis[j,i] = dis[i,j]}
            if(method == 'ward' ){
              mss = rbind(colMeans(clu_mtx[1:a,,i]), colMeans(clu_mtx[1:b,,j]))
              dis[i,j] = as.matrix(dist(mss))[1,2]
              dis[i,j] = dis[i,j] * sqrt(a*b/(a+b))
              dis[j,i] = dis[i,j]
            }
            if(method == 'tv' ){
              dis[i,j] = tvd(colMeans(clu_mtx[1:a,,i]), colMeans(clu_mtx[1:b,,i]))
              dis[j,i] = dis[i,j]
            }
            if(method == 'cid' ){
              dis[i,j] = diss.CID(colMeans(clu_mtx[1:a,,i]), colMeans(clu_mtx[1:b,,i]))
              dis[j,i] = dis[i,j]
            }
            
          } 
        }
      }
    }
  }
  return(list(n_element_in_each_clu=n_element_in_each_clu, res=res, min_dis=min_dis))
}

# EACH SAMPLE IS AN INITIAL CLUSTER
robustcluster_1 = function(x, n, method,tau){
  dimension = dim(x)
  n_rep = dimension[1]
  n_freq = dimension[2]
  n_clu = dimension[3]
  min_dis = rep(0, n_clu)
  res = 1:n_clu
  n_element_in_each_clu = matrix(1, 1, n_clu)
  clu_mtx = array(0,c(n_clu*n_rep, n_freq,n_clu))
  for (i in 1:n_clu){
    clu_mtx[1:n_rep, , i] = x[, , i]
  }
  if(method == 'cr'){dis = distcr(x,tau)}
  if(method == 'ms'){dis = distcr2(x,tau)}
  
  if(method == 'tv'){
    if(n_rep>1){
      dis=matrix(0,n_clu,n_clu)
      for(i in 1:(n_clu-1)){
        for(j in (i+1):n_clu){
          dis[i,j]=tvd(colMeans( x[,,i]),colMeans( x[,,j]))
          dis[j,i]=dis[i,j]
        }
      }}
    if(n_rep==1){
      dis=matrix(0,n_clu,n_clu)
      for(i in 1:(n_clu-1)){
        for(j in (i+1):n_clu){
          dis[i,j]=tvd(x[1,,i],x[1,,j])
          dis[j,i]=dis[i,j]
        }}
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  if(method == 'fm'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:(n_clu-1)){
      for(j in (i+1):n_clu){
        dis[i,j]=diss(fmed(x[,,i]),fmed(x[,,j]))
        dis[j,i]=dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  if(method == 'mean'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:n_clu){
      for(j in 1:n_clu){
        ms = rbind(colMeans(x[,,i]), colMeans(x[,,j]))
        dis[i,j] = as.matrix(dist(ms))[1,2]
        dis[j,i]=dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  if(method == 'ward'){
    dis=matrix(0,n_clu,n_clu)
    for(i in 1:n_clu){
      for(j in 1:n_clu){
        ms = rbind(colMeans(x[,,i]), colMeans(x[,,j]))
        dis[i,j] = as.matrix(dist(ms))[1,2] * sqrt(n_rep^2/(2*n_rep))
        dis[j,i] = dis[i,j]
      }
    }
    for(i in 1:n_clu){
      dis[i,i] = 10000000
    }
  }
  
  
  rcd = 1
  for(k in n_clu:(n+1)){
    index_min = which(dis==min(dis), arr.ind=TRUE)
    index_min = index_min[1,]
    if(index_min[2] <= index_min[1]){
      temp = index_min[1]
      index_min[1] = index_min[2]
      index_min[2] = temp}
    min_dis[rcd] = dis[index_min[1], index_min[2]]
    rcd = rcd + 1
    clu_mtx[(n_rep*n_element_in_each_clu[index_min[1]]+1):(n_rep*(n_element_in_each_clu[index_min[1]]+n_element_in_each_clu[index_min[2]])), , index_min[1]]=clu_mtx[1:(n_rep*n_element_in_each_clu[index_min[2]]), , index_min[2]]
    clu_mtx[ , ,index_min[2]] = matrix(0,n_clu*n_rep, n_freq)
    n_element_in_each_clu[index_min[1]] = n_element_in_each_clu[index_min[1]]+n_element_in_each_clu[index_min[2]]
    n_element_in_each_clu[index_min[2]] = 0
    index = which(res==res[index_min[2]])
    res[index] = res[index_min[1]]
    
    for(i in 1:(n_clu-1)){
      for (j in (i+1):n_clu){
        a = notzero(clu_mtx[, , i])
        b = notzero(clu_mtx[, , j])
        if (a==0 || b==0){
          dis[i,j] = 10000000
          dis[j,i] = 10000000}
        else{
          if(i==index_min[1]){
            if(method == 'cr'){
              if(a >=7 & b >=7){
                dis[i,j] = ddcr(clu_mtx[1:a, , i], clu_mtx[1:b, , j], tau)
                dis[j,i] = dis[i,j]}
              else{
                mss = rbind(colMeans(matrix(clu_mtx[1:a,,i],a,n_freq)), colMeans(matrix(clu_mtx[1:b,,i],b,n_freq)))
                dis[i,j] = as.matrix(dist(mss))[1,2]
                dis[i,j] = dis[i,j] * sqrt(a*b/(a+b))
                dis[j,i] = dis[i,j]}}
            
            if(method == 'ms'){
              if(a >=15 & b >=15){
                dis[i,j] = ddcr2(clu_mtx[1:a, , i], clu_mtx[1:b, , j], tau)
                dis[i,j] = dis[i,j] * sqrt(a*b/(a+b))
                dis[j,i] = dis[i,j]}
              else{
                mss = rbind(colMeans(matrix(clu_mtx[1:a,,i],a,n_freq)), colMeans(matrix(clu_mtx[1:b,,i],b,n_freq)))
                dis[i,j] = as.matrix(dist(mss))[1,2]
                dis[i,j] = dis[i,j] * sqrt(a*b/(a+b))
                dis[j,i] = dis[i,j]}}
            if(method == 'fm'){
              dis[i,j]=diss(fmed(clu_mtx[1:a,,i]), fmed(clu_mtx[1:b,,j]))
              dis[j,i] = dis[i,j]}
            if(method == 'mean'){
              mss = rbind(colMeans(matrix(clu_mtx[1:a,,i],a,n_freq)), colMeans(matrix(clu_mtx[1:b,,i],b,n_freq)))
              dis[i,j] = as.matrix(dist(mss))[1,2]
              dis[j,i] = dis[i,j]}
            if(method == 'ward' ){
              mss = rbind(colMeans(matrix(clu_mtx[1:a,,i],a,n_freq)), colMeans(matrix(clu_mtx[1:b,,i],b,n_freq)))
              dis[i,j] = as.matrix(dist(mss))[1,2]
              dis[i,j] = dis[i,j] * sqrt(a*b/(a+b))
              dis[j,i] = dis[i,j]
            }
            if(method == 'tv' ){
              dis[i,j] = tvd(colMeans(matrix(clu_mtx[1:a,,i],a,n_freq)), colMeans(matrix(clu_mtx[1:b,,i],b,n_freq)))
              dis[j,i] = dis[i,j]
            }
          } 
        }
      }
    }
  }
  return(list(n_element_in_each_clu=n_element_in_each_clu, res=res, min_dis=min_dis))
}


#standerize the clusters
stdclu=function(x){
  len = length(x)
  s = rep(0,len)
  cl = rep(0,len)
  for(i in 1:len){
    cl[i] = (x[i] == i)
  }
  num = sum(cl)
  t = which(cl == 1)
  for(i in 1:num){
    p = which(x == t[i])
    s[p] = i
  }
  return(s)
}

####rank sum test
ranktest = function(dataX, dataY, dataRef)
{
  n = dim(dataX)[2];
  m = dim(dataY)[2];
  r = dim(dataRef)[2];
  order = integer(n+m);
  for(i in 1:m)
  {
    sample = cbind(dataRef, dataY[,i]);
    result = fbplot(sample, plot=F);
    order[i] = sum(result$depth[1:r] <= result$depth[r+1])
  }
  for(i in 1:n)
  {
    sample = cbind(dataRef,dataX[,i]);
    result = fbplot(sample, plot=F);
    order[i+m] = sum(result$depth[1:r] <= result$depth[r+1])
  }
  rk = sort.int(order, index.return = T)$ix;
  w = sum(rk[1:m]);
  return(w)
}

#combination
combinat = function(n,p){
  if (n<p){combinat=0}
  else {combinat = exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}

# MBD
fBD2 = function(data){
  p = dim(data)[1]
  n = dim(data)[2]
  rmat = apply(data, 1, rank)
  down = apply(rmat, 1, min)-1
  up = n-apply(rmat, 1, max)
  (up*down+n-1) / combinat(n,2)
}
########MBD########
fMBD = function(data){
  p = dim(data)[1]
  n = dim(data)[2]
  rmat = apply(data,1,rank)
  down = rmat-1
  up = n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}

#L2 DIST
diss=function(x1, x2){ return(sqrt(sum((x1 - x2) ^ 2))) }

# CONTOUR 2D
c2d = function(xdata, ydata, tau, xli=0, yli=0 ,grid_size=1, plot=T){
  t_length = 1
  t = seq(0.1,1,length.out = t_length)
  n_each = length(xdata)
  Y1=matrix(0, nrow=length(t), ncol = n_each); 
  Y2=matrix(0, nrow=length(t), ncol = n_each); 
  Y1[1, ] = 1 + xdata
  Y2[1, ] = 0.5 + ydata
  m = grid_size^2
  Y1_alltime = vec(t(Y1)) 
  Y2_alltime = vec(t(Y2)) 
  Y = cbind(Y1_alltime, Y2_alltime)
  n = n_each*length(t)
  df = 7
  basis_model = bs(x = t, df= df, degree = 3,
                   Boundary.knots = c(0,1), intercept = TRUE)
  vec_n= vector()
  for(i in 1:t_length){
    vec_i = c(rep(basis_model[i,], times= n_each))
    vec_n = c(vec_n, vec_i)
  }
  X =  matrix(vec_n, nrow = n, ncol =df, byrow=T)
  nu = rep(1/n, n)
  mu = rep(1/m, m)
  x <- seq(0, 1, length.out = grid_size)
  y <- seq(0, 1, length.out = grid_size)
  U <- as.matrix(expand.grid(x = x, y = y))
  obj1 = vec(as.matrix(nu))
  obj2 = vec(mu%*%t(nu)%*%X)
  model_dual = list()
  
  model_dual$obj = c(obj1, obj2)
  
  A1_dual = Matrix(kronecker(X= diag(n), Y= rep(1,m)), sparse = T)
  A2_dual = Matrix(kronecker(X= X , Y= diag(m)), sparse = T)
  model_dual$A = cbind(A1_dual,A2_dual)
  
  rhs_dual = vec(U%*%t(Y))
  model_dual$rhs <- rhs_dual
  model_dual$vtype <- 'C'  
  model_dual$sense <- '>'
  model_dual$modelsense <- 'min'
  
  params <- list(OutputFlag=0)
  result_dual <- gurobi(model_dual, params)
  psi = result_dual$x[1:n]
  vec_b = result_dual$x[-(1:n)]
  b = matrix(vec_b, nrow = m, ncol = df) 
  
  library(numDeriv)
  covariate = predict(basis_model, newx = t[1])[1,]
  bTx= apply(b, 1, function(x) t(x)%*%covariate)
  
  finite.differences <- function(x, y) {
    
    stepsize = abs(x[1,1]-x[2,1])
    n <- dim(x)[1]
    
    # Initialize a vector of length n to enter the derivative approximations
    fdx_u <- vector(length = grid_size)
    fdx_f = rep(0,grid_size); fdx_b = rep(0,grid_size)
    # Iterate through the values using the forward differencing method
    fdx = matrix(0, nrow=grid_size, ncol=grid_size)
    
    for(u in 1:grid_size){
      for (i in 1:(grid_size-1)) {
        fdx_f[i] <- (y[i+1 +grid_size*(u-1)] - y[i +grid_size*(u-1)]) / stepsize
      }
      
      # Iterate through the values using the backward differencing method
      for (i in 2:grid_size) {
        fdx_b[i] <- (y[i +grid_size*(u-1)] - y[i-1 +grid_size*(u-1)]) / stepsize
      }
      
      fdx_u <- (fdx_f+fdx_b)/2
      fdx_u[1] <- fdx_f[1]
      fdx_u[grid_size] <- fdx_b[grid_size]
      
      fdx[,u] = fdx_u
    }
    return(fdx)
  }
  
  beta1Tx = finite.differences(U, bTx)
  fmatrix= matrix(bTx, nrow=grid_size, ncol=grid_size)
  b2Tx = vec(t(fmatrix))
  
  beta2Tx = t(finite.differences(U, b2Tx))
  Y_cap =cbind(c(beta1Tx), c(beta2Tx))
  norm_vec <- function(x) sqrt(sum(x^2))
  norm_max = function(x) max(abs(x[1]), abs(x[2]))
  U_centered = U-0.5
  radius = seq(0.00001, 1, length.out = grid_size)
  
  num <- grid_size # number of points you want on the unit circle
  pts_all <- t(sapply(1:num,function(p) c(radius[1]*cos(2*p*pi/num),radius[1]*sin(2*p*pi/num))))
  for(i in 2:grid_size){
    pts.circle <- t(sapply(1:num,function(p) c(radius[i]*cos(2*p*pi/num),radius[i]*sin(2*p*pi/num))))
    pts_all = rbind(pts_all, pts.circle)
  }
  U_n = pts_all
  
  library(pdist)
  dist_mat = as.matrix(pdist(X=U_n, Y=U_centered))
  ranks = rank(c(dist_mat))
  cmat = matrix(ranks, nrow=grid_size^2, ncol=grid_size^2)
  library(adagio)
  matching = assignment(cmat)
  perm = matching$perm   
  
  library(alphahull)
  U_number = which(apply(U_n, 1, norm_vec) < tau)
  U_centered_corresponding = perm[U_number]
  Y_tau = Y_cap[U_centered_corresponding,]
  ahull.obj = ahull(Y_tau[,1]-1, Y_tau[,2]-0.5, 800)
  edges = ahull.obj$ashape.obj$edges
  ind = ahull.obj$ashape.obj$alpha.extremes
  pol = point.in.polygon(xdata+1,ydata+0.5,Y_tau[ind,1],Y_tau[ind,2])
  if(plot == T){
    if(length(xli) == 1){
      plot(Y1[1,]-1, Y2[1,]-0.5, xlab='MO', ylab='VO', pch=20, col='gray',xaxt='n',yaxt='n')
      plot(ahull.obj,lwd=c(1,1,2), add=T, wpoints=F)
      points(xdata[which(pol==1)], ydata[which(pol==1)],col='gray31', pch=20)}
    
    if(length(xli) == 2){
      plot(Y1[1,]-1, Y2[1,]-0.5, xlab='MO', ylab='VO', pch=20, col='gray', xlim=xli, ylim=yli,xaxt='n',yaxt='n')
      plot(ahull.obj,lwd=c(1,1,2), add=T, wpoints=F)
      points(xdata[which(pol==1)], ydata[which(pol==1)],col='gray31',pch=20)}
  }
  
  return(list(edges_index=ind, pts_in_pol=which(pol==1)))
}

#MS
hardin_factor_numeric <- function(n, dimension){
  h <- floor((n+dimension+1)/2)
  alpha <- (n-h)/n
  q_alpha <- qchisq(1-alpha, dimension)
  c_alpha <- (1 - alpha)/pchisq(q_alpha, dimension + 2)
  c2 <- -pchisq(q_alpha, dimension+2)/2
  c3 <- -pchisq(q_alpha, dimension + 4)/2
  c4 <- 3*c3
  b1 <- c_alpha*(c3-c4)/(1-alpha)
  b2 <- 0.5 + c_alpha/(1-alpha)*(c3-q_alpha*(c2+(1-alpha)/2)/dimension)
  v1  <- (1-alpha)*(b1^2)*
    (alpha*(c_alpha*q_alpha/dimension-1)^2-1)-2*c3*c_alpha^2*
    (3*(b1-dimension*b2)^2+(dimension+2)*b2*(2*b1-dimension*b2))
  v2 <- n*(b1*(b1-dimension*b2)*(1-alpha))^2*c_alpha^2
  v <- v1/v2
  m_asy <- 2/(c_alpha^2*v)
  m <- m_asy*exp(0.725-0.00663*dimension-0.078*log(n))
  if (m < dimension){ #if m is >= dimension, then line 18 works, if not change m to m_asy
    m <- m_asy
  }
  a1 <- rchisq(10000,dimension + 2)
  a2 <- rchisq(10000,dimension, h/n)
  c <- sum(a1 < a2)/(10000*h/n)
  factors <- c * (m - dimension + 1)/(dimension * m)
  cutoff <- qf(0.993, dimension, m - dimension + 1)
  list(factor1 = factors, factor2 = cutoff)
}

msplot1=function (dts, data_depth = c("random_projections"), n_projections = 200, 
                  seed = NULL, return_mvdir = TRUE, plot = TRUE, plot_title = "Magnitude Shape Plot", 
                  title_cex = 1.5, show_legend = T, ylabel = "VO", xlabel) {
  data_dim <- dim(dts)
  n <- data_dim[1]
  dir_result <- dir_out(dts, data_depth = data_depth, n_projections = n_projections, 
                        seed = seed)
  if (length(data_dim) == 2) {
    dist <- dir_result$distance
    rocke_factors <- hardin_factor_numeric(n, 2)
    rocke_factor1 <- rocke_factors$factor1
    rocke_cutoff <- rocke_factors$factor2
    cutoff_value <- rocke_cutoff/rocke_factor1
    outliers_index <- which(dist > cutoff_value*2.5)
    median_curve <- which.min(dist)
    if (plot) {
      myx <- dir_result$mean_outlyingness
      myy <- dir_result$var_outlyingness
      if (missing(xlabel)) 
        xlabel <- "MO"
    }
  }
  else if (length(data_dim) == 3) {
    d <- data_dim[3]
    rocke_factors <- hardin_factor_numeric(n = n, dimension = d + 
                                             1)
    rocke_factor1 <- rocke_factors$factor1
    rocke_cutoff <- rocke_factors$factor2
    cutoff_value <- rocke_cutoff/rocke_factor1
    outliers_index <- which(dir_result$distance > (cutoff_value))
    median_curve <- which.min(dir_result$distance)
    if (plot) {
      myx <- sqrt(rowSums(dir_result$mean_outlyingness^2, 
                          na.rm = T))
      myy <- dir_result$var_outlyingness
      if (missing(xlabel)) 
        xlabel <- "||MO||"
    }
  }
  if (plot) {
    plot(myx, myy, type = "n", xlab = xlabel, ylab = ylabel, 
         xlim = range(myx) + c(-sd(myx), 1.5 * sd(myx)), ylim = range(myy) + 
           c(-0.2 * sd(myy), 1 * sd(myy)), axes = F, col.lab = "gray20")
    axis(1, col = "white", col.ticks = "grey61", lwd.ticks = 0.5, 
         tck = -0.025, cex.axis = 0.9, col.axis = "gray30")
    axis(2, col = "white", col.ticks = "grey61", lwd.ticks = 0.5, 
         tck = -0.025, cex.axis = 0.9, col.axis = "gray30")
    grid(col = "grey75", lwd = 0.3)
    box(col = "grey51")
    if (length(outliers_index > 0)) {
      points(myx[-outliers_index], myy[-outliers_index], 
             bg = "gray60", pch = 21)
      points(myx[outliers_index], myy[outliers_index], 
             pch = 3)
    }
    else {
      points(myx, myy, bg = "gray60", pch = 21)
    }
    mtext(plot_title, 3, adj = 0.5, line = 1, cex = title_cex, 
          col = "gray20")
    if (show_legend) {
      legend("topright", legend = c("normal", "outlier"), 
             pch = c(21, 3), cex = 1, pt.bg = "gray60", col = "gray0", 
             text.col = "gray30", bty = "n", box.lwd = 0.1, 
             xjust = 0, inset = 0.01)
    }
  }
  if (return_mvdir) {
    return(list(outliers = outliers_index, median_curve = median_curve, 
                mean_outlyingness = dir_result$mean_outlyingness, 
                var_outlyingness = dir_result$var_outlyingness))
  }
  else {
    return(list(outliers = outliers_index, median_curve = median_curve))
  }
}

#PLOT A DISTANCE MATRIX
plotdis = function(dis, name=''){
  for(i in 1:dim(dis)[1]){
    for(j in 1:dim(dis)[2]){
      if(dis[i,j]==10000000){
        dis[i,j]=0
      }
    }
  }
  dis = dis/max(dis)
  image.plot(1:dim(dis)[1],1:dim(dis)[1],dis,xlab='',ylab='', main=name )
  axis(side = 1, tck = 0.00) ;axis(side = 2, tck = 0.00)
  return(dis)
}

#PLOT A CLUSTER
plotclu = function(x, y=0, z=0){
  lowx=min(x)
  hix=max(x)
  if(length(dim(y))!=0){
    lowy=min(y)
    hiy=max(y)
  }
  else{
    lowy=lowx+1
    hiy=hix-1
  }
  
  if(length(dim(z))!=0){
    lowz = min(z)
    hiz = max(z)
  }
  else{
    lowz = lowx+1
    hiz = hix-1
  }
  low = min(lowx, lowy, lowz)
  hi = max(hix, hiy, hiz)
  
  
  n = dim(x)[1]
  plot(x[1,],type='l', col=1, ylim=c(low, hi))
  for(i in 2:n){
    lines(x[i,],col=1)
  }
  if(length(dim(y))!=0){
    for(i in 1:dim(y)[1]){
      lines(y[i,],col=2)
    }
  }
  if(length(dim(z))!=0){
    for(i in 1:dim(z)[1]){
      lines(z[i,],col=3)
    }
  }
}

#GENERATE GAUSSIAN PROCESS
gaussprocess <- function(n, K = function(s, t) {min(s, t)}) {
t <- seq(from = 0, to = 1, length.out = n)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0,  n), Sigma = Sigma)
  path <- path   # Must always start at "start"
  
  return(path)
}

# SGN
sgn = function(pro){
  if(pro>=0.5) return(1)
  if(pro<0.5)  return(-1)
}

# AR SPECTRUM
specar = function(ar, l){
  p = length(ar)
  omega = seq(from=-0.5, to=0.5, length.out=l+3)
  omega = omega[2:(length(omega)-1)]
  sp = matrix(0,p,l+1)
  for(j in 1:p){
    sp[j, ] = ar[j] * exp( -2 * (j) * pi * (0+1i) * omega)
  }
  sp = colSums(sp)
  spe = 1 / (1 - sp)
  return((abs(spe)) ^ 2)
}


LLR = function(x, y){
  ratio = exp(x)/exp(y)
  return(sum(0.5 * log(ratio) + 0.5 * log(1 / ratio)))
}

mi = function(x,y){
  xymin = rep(0,length(x))
  for(i in 1:length(x)){
    xymin[i] = min(x[i],y[i])
  }
  return(xymin)
}


flat = function(x){
  p = dim(x)[1]
  q = dim(x)[2]
  r = dim(x)[3]
  mat = matrix(0, p*r, q)
  for(i in 1:r){
    mat[((i-1)*p+1):(i*p),] = x[,,i]}
  return(mat)
}


flatm = function(x){
  p = dim(x)[1]
  q = dim(x)[2]
  r = dim(x)[3]
  mat = matrix(0, r, q)
  for(i in 1:r){
    mat[i,] = colMeans(x[,,i])}
  return(mat)
}



tclust1=function(x, nc, rate){
  x = flat(x)
  res = tclust(x,nc,alpha=rate)
  return(res$cluster)
}


