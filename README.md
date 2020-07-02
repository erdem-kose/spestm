# SPeST: Spectral Estimation Library


In readme, for usage example I will inspect the popular spectral estimation methods with audio file consists G5 note with trumpet. Each of models will be shown in whole and short-time structure. My personal favorite is Min-Norm method, because it gives Autoregressive parameters and a good Power Spectral Density at the same time.


# Data Load

```matlab:Code
clear; clc; close all;

[x,fs]=audioread('trumpet-G5.wav');
N=size(x,1);
d=size(x,2);
t=fs*((1:N)-1)';

xzm=x-mean(x,1);

x_channel_labels=cell(d,1);
for i=1:d
    x_channel_labels{i}=['Channel ' num2str(i)];
end
```

### Short-time Partitioning

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;N_{\textrm{partition}}&space;=\textrm{ceil}\left(\frac{N}{\textrm{Partition}\;\textrm{Count}}\right),N_{\textrm{new}}&space;=\textrm{Partition}\;\textrm{Count}*N_{\textrm{partition}}&space;\\&space;x_{\textrm{new}}&space;=\left\lbrack&space;\begin{array}{cc}&space;x&space;&&space;\textrm{zerox}\left(N_{\textrm{new}}&space;-N,d\right)&space;\end{array}\right\rbrack&space;\\&space;x_{\textrm{shorttime}}&space;=\textrm{reshape}\left(x,\left\lbrack&space;N_{\textrm{partition}}&space;,\textrm{Partition}\;\textrm{Count}\right\rbrack&space;\right)&space;\end{array}"/>


```matlab:Code
partition_count=500;
partition_N=ceil(N/partition_count);
Nnew=partition_count*partition_N;
x_st=[x; zeros(Nnew-N,d)];
x_st=reshape(x_st,[partition_N,partition_count]);
x_st=x_st-mean(x_st,1);

x_stax = [0 Nnew/fs]; %t axis
y_stax = [0 fs/2]; %f axis
```

### Plot Time Series

```matlab:Code
plot(t,x)
ylabel('Amplitude ');
xlabel('Time(Sec)');
legend(x_channel_labels)
grid minor;
```


![figure_0.png](readme_images/figure_0.png)

# Spectral Estimation Library

```matlab:Code
spestm_obj=spestm_lib();
spestm_obj.f_range='half';
spestm_obj.window_type='hann';
spestm_obj.x=xzm;
%spestm_obj.x=spestm_obj.side_awgn(20);
spestm_obj.fs=fs;
spestm_obj.p=24;
spestm_obj.q=4;
spestm_obj.M=96;
spestm_obj=spestm_obj.init();
```


```matlab:Code
spestm_obj_st=spestm_lib();
spestm_obj_st.f_range='half';
spestm_obj_st.window_type='hann';
spestm_obj_st.x=x_st;
spestm_obj_st.x=spestm_obj_st.side_awgn(20);
spestm_obj_st.fs=fs;
spestm_obj_st.p=12;
spestm_obj_st.q=2;
spestm_obj_st.M=24;
spestm_obj_st=spestm_obj_st.init();
```

## Periodogram PSD


The most basic one is Periodogram, it's just amplitude spectrum of autocorrelation function <img src="https://latex.codecogs.com/gif.latex?\inline&space;r_X&space;\left\lbrack&space;k\right\rbrack"/>'s Fouirer transform.



<img src="https://latex.codecogs.com/gif.latex?x_w&space;\left\lbrack&space;n\right\rbrack&space;=w\left\lbrack&space;n\right\rbrack&space;x\left\lbrack&space;n\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?r_X&space;\left\lbrack&space;k\right\rbrack&space;=\frac{1}{f_{s\;}&space;}\sum_{n=0}^{N-1}&space;x_w&space;\left\lbrack&space;n\right\rbrack&space;\bar{x_w&space;}&space;\left\lbrack&space;n-k\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?P_X^{\textrm{per}}&space;\left(f\right)=\left|\;\sum_{k=0}^{N-1}&space;r_X&space;\left\lbrack&space;k\right\rbrack&space;e^{-\textrm{j2}\pi&space;\frac{f}{f_s&space;}n}&space;\right|"/>


```matlab:Code
[y_per,f]=spestm_obj.psd_periodogram();
y_per_st=spestm_obj_st.psd_periodogram();
```


```matlab:Code
semilogy(f,y_per);
ylabel('PSD(W/Hz^2) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_1.png](readme_images/figure_1.png)


```matlab:Code
image_y=log10(y_per_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_2.png](readme_images/figure_2.png)

## Blackman-Tukey PSD


It's the variant of Periodogram, but differs in windowing. We window <img src="https://latex.codecogs.com/gif.latex?\inline&space;R_X&space;\left\lbrack&space;k\right\rbrack"/> instead of <img src="https://latex.codecogs.com/gif.latex?\inline&space;x\left\lbrack&space;n\right\rbrack"/>.



<img src="https://latex.codecogs.com/gif.latex?r_X&space;\left\lbrack&space;k\right\rbrack&space;=\frac{1}{f_{s\;}&space;}\sum_{n=0}^{N-1}&space;x\left\lbrack&space;n\right\rbrack&space;\bar{x}&space;\left\lbrack&space;n-k\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?r_{X_w&space;}&space;\left\lbrack&space;k\right\rbrack&space;=w\left\lbrack&space;k\right\rbrack&space;R_X&space;\left\lbrack&space;k\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?P_X^{\textrm{bt}}&space;\left(f\right)=\left|\;\sum_{k=0}^{N-1}&space;r_{X_w&space;}&space;\left\lbrack&space;k\right\rbrack&space;e^{-\textrm{j2}\pi&space;\frac{f}{f_s&space;}n}&space;\right|"/>


```matlab:Code
[y_bt,f]=spestm_obj.psd_blackmantukey();
y_bt_st=spestm_obj_st.psd_blackmantukey();
```


```matlab:Code
semilogy(f,y_bt);
ylabel('PSD(W/Hz) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_3.png](readme_images/figure_3.png)


```matlab:Code
image_y=log10(y_bt_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_4.png](readme_images/figure_4.png)

## Capon PSD


From now on, we don't mention about windowing. We will use input signal directly. Capon PSD, is one of the best because of detail parameter <img src="https://latex.codecogs.com/gif.latex?\inline&space;g"/>. If <img src="https://latex.codecogs.com/gif.latex?\inline&space;g"/> increases, details of noisy input signal can be observed clearly; but it has a limit naturally.



<img src="https://latex.codecogs.com/gif.latex?r_X&space;\left\lbrack&space;k\right\rbrack&space;=\frac{1}{f_{s\;}&space;}\sum_{n=0}^{N-1}&space;x\left\lbrack&space;n\right\rbrack&space;\bar{x}&space;\left\lbrack&space;n-k\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?R_X&space;=\left\lbrack&space;\begin{array}{cccc}&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-1\right)\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-2\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;p-1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;p-2\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;\end{array}\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?{\mathit{\mathbf{e}}}_{\textrm{vect}}&space;={\left\lbrack&space;\begin{array}{cccc}&space;1&space;&&space;e^{-\textrm{j2}\pi&space;f\;}&space;&space;&&space;\cdots&space;&space;&&space;e^{-\textrm{j2}\pi&space;f\left(p-1\right)\;}&space;&space;\end{array}\right\rbrack&space;}^T"/>


<img src="https://latex.codecogs.com/gif.latex?P_X^{\textrm{capon}}&space;\left(f\right)=\left|\frac{{\mathit{\mathbf{e}}}_{\textrm{vect}}^{\mathit{\mathbf{T}}}&space;\left(R_X^{-\left(g-1\right)}&space;\right){\mathit{\mathbf{e}}}_{\textrm{vect}}&space;}{\;{\mathit{\mathbf{e}}}_{\textrm{vect}}^{\mathit{\mathbf{T}}}&space;\left(R_X^{-\left(g\right)}&space;\right){\mathit{\mathbf{e}}}_{\textrm{vect}}&space;}\right|"/>


```matlab:Code
[y_cap,f]=spestm_obj.psd_capon();
y_cap_st=spestm_obj_st.psd_capon();
```


```matlab:Code
semilogy(f,y_cap);
ylabel('PSD(W/Hz^2) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_5.png](readme_images/figure_5.png)


```matlab:Code
image_y=log10(y_cap_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_6.png](readme_images/figure_6.png)

## Autoregressive (Yule-Walker) PSD

<img src="https://latex.codecogs.com/gif.latex?r_X&space;\left\lbrack&space;k\right\rbrack&space;=\frac{1}{f_{s\;}&space;}\sum_{n=0}^{N-1}&space;x\left\lbrack&space;n\right\rbrack&space;\bar{x}&space;\left\lbrack&space;n-k\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?R_X&space;=\left\lbrack&space;\begin{array}{cccc}&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-1\right)\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-2\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;p-1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;p-2\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;\end{array}\right\rbrack"/>



Solution of the equation above will give parameter vector <img src="https://latex.codecogs.com/gif.latex?\inline&space;a"/>.



<img src="https://latex.codecogs.com/gif.latex?\left\lbrack&space;\begin{array}{cccc}&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-1\right)\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-2\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;p-1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;p-2\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;\end{array}\right\rbrack&space;\left\lbrack&space;\begin{array}{c}&space;a\left\lbrack&space;1\right\rbrack&space;\\&space;a\left\lbrack&space;2\right\rbrack&space;\\&space;\vdots&space;\\&space;a\left\lbrack&space;p\right\rbrack&space;&space;\end{array}\right\rbrack&space;=-\left\lbrack&space;\begin{array}{c}&space;r_X&space;\left\lbrack&space;1\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;2\right\rbrack&space;\\&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;p\right\rbrack&space;&space;\end{array}\right\rbrack"/>



Parameter <img src="https://latex.codecogs.com/gif.latex?\inline&space;\sigma&space;{\;}^2"/> can be obtained by after obtaining parameter vector <img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathit{\mathbf{a}}={\left\lbrack&space;\begin{array}{cccc}
1&space;&&space;a\left\lbrack&space;1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;a\left\lbrack&space;p\right\rbrack&space;
\end{array}\right\rbrack&space;}^T"/>



<img src="https://latex.codecogs.com/gif.latex?\left\lbrack&space;\begin{array}{cccc}&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p\right)\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;-\left(p-1\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;p\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;p-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;0\right\rbrack&space;&space;\end{array}\right\rbrack&space;\left\lbrack&space;\begin{array}{c}&space;1\\&space;a\left\lbrack&space;1\right\rbrack&space;\\&space;\vdots&space;\\&space;a\left\lbrack&space;p\right\rbrack&space;&space;\end{array}\right\rbrack&space;=-\left\lbrack&space;\begin{array}{c}&space;\sigma&space;{\;}^2&space;\\&space;0\\&space;\vdots&space;\\&space;0&space;\end{array}\right\rbrack"/>



Then using these parameters



<img src="https://latex.codecogs.com/gif.latex?P_X^{\textrm{ARyw}}&space;\left(f\right)=\frac{\sigma&space;{\;}^2&space;}{\;\left|{\mathit{\mathbf{e}}}_{\textrm{vect}}^{\mathit{\mathbf{T}}}&space;\mathit{\mathbf{a}}{\left|\right.}^2&space;\right.}"/>


<img src="https://latex.codecogs.com/gif.latex?{\mathit{\mathbf{e}}}_{\textrm{vect}}&space;={\left\lbrack&space;\begin{array}{cccc}&space;1&space;&&space;e^{-\textrm{j2}\pi&space;f\;}&space;&space;&&space;\cdots&space;&space;&&space;e^{-\textrm{j2}\pi&space;\textrm{fp}\;}&space;&space;\end{array}\right\rbrack&space;}^T"/> 


```matlab:Code
[y_aryw,f] = spestm_obj.psd_aryulewalker();
y_aryw_st=spestm_obj_st.psd_aryulewalker();
```


```matlab:Code
semilogy(f,y_aryw);
ylabel('PSD(W/Hz) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_7.png](readme_images/figure_7.png)


```matlab:Code
image_y=log10(y_aryw_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_8.png](readme_images/figure_8.png)

## Autoregressive (Modified Covariance) PSD

```matlab:Code
[y_armc,f] = spestm_obj.psd_armodcov();
y_armc_st=spestm_obj_st.psd_armodcov();
```


```text:Output
Rank deficient, p is changed to 0
Rank deficient, p is changed to 0
Rank deficient, p is changed to 0
```


```matlab:Code
semilogy(f,y_armc);
ylabel('PSD(W/Hz) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_9.png](readme_images/figure_9.png)


```matlab:Code
image_y=log10(y_armc_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_10.png](readme_images/figure_10.png)

## Autoregressive (Burg) PSD

```matlab:Code
[y_arburg,f]=spestm_obj.psd_arburg();
y_arburg_st=spestm_obj_st.psd_arburg();
```


```matlab:Code
semilogy(f,y_arburg);
ylabel('PSD(W/Hz^2) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_11.png](readme_images/figure_11.png)


```matlab:Code
image_y=log10(y_arburg_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_12.png](readme_images/figure_12.png)

## Autoregressive Moving Average PSD

<img src="https://latex.codecogs.com/gif.latex?r_X&space;\left\lbrack&space;k\right\rbrack&space;=\frac{1}{f_{s\;}&space;}\sum_{n=0}^{N-1}&space;x\left\lbrack&space;n\right\rbrack&space;\bar{x}&space;\left\lbrack&space;n-k\right\rbrack"/>


<img src="https://latex.codecogs.com/gif.latex?R_X^{\textrm{modified}}&space;=\left\lbrack&space;\begin{array}{cccc}&space;r_X&space;\left\lbrack&space;q+1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;q\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;q-\left(p-2\right)\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;q+2\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;q+1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;q-\left(p-3\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;q+p\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;q+p-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;q+1\right\rbrack&space;&space;\end{array}\right\rbrack"/>



Solution of the equation above will give parameter vector <img src="https://latex.codecogs.com/gif.latex?\inline&space;a"/>.



<img src="https://latex.codecogs.com/gif.latex?\left\lbrack&space;\begin{array}{cccc}&space;r_X&space;\left\lbrack&space;q+1\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;q\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;q-\left(p-2\right)\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;q+2\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;q+1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;q-\left(p-3\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;q+p\right\rbrack&space;&space;&&space;r_X&space;\left\lbrack&space;q+p-1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;r_X&space;\left\lbrack&space;q+1\right\rbrack&space;&space;\end{array}\right\rbrack&space;\left\lbrack&space;\begin{array}{c}&space;a\left\lbrack&space;1\right\rbrack&space;\\&space;a\left\lbrack&space;2\right\rbrack&space;\\&space;\vdots&space;\\&space;a\left\lbrack&space;p\right\rbrack&space;&space;\end{array}\right\rbrack&space;=-\left\lbrack&space;\begin{array}{c}&space;r_X&space;\left\lbrack&space;q+2\right\rbrack&space;\\&space;r_X&space;\left\lbrack&space;q+3\right\rbrack&space;\\&space;\vdots&space;\\&space;r_X&space;\left\lbrack&space;q+p+1\right\rbrack&space;&space;\end{array}\right\rbrack"/>



Filtering input signal <img src="https://latex.codecogs.com/gif.latex?\inline&space;x"/> with parameter vector <img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathit{\mathbf{a}}={\left\lbrack&space;\begin{array}{cccc}
1&space;&&space;a\left\lbrack&space;1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;a\left\lbrack&space;p\right\rbrack&space;
\end{array}\right\rbrack&space;}^T"/>



<img src="https://latex.codecogs.com/gif.latex?x_f&space;\left\lbrack&space;n\right\rbrack&space;=\sum_{k=0}^p&space;x\left\lbrack&space;n\right\rbrack&space;a\left\lbrack&space;k-n\right\rbrack"/>



Find parameter vector <img src="https://latex.codecogs.com/gif.latex?\inline&space;{\mathit{\mathbf{a}}}_{\mathit{\mathbf{f}}}&space;={\left\lbrack&space;\begin{array}{cccc}
1&space;&&space;a_f&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;a_f&space;\left\lbrack&space;p\right\rbrack&space;
\end{array}\right\rbrack&space;}^T"/> for  <img src="https://latex.codecogs.com/gif.latex?\inline&space;x_f&space;\left\lbrack&space;n\right\rbrack"/> with Yule-Walker equations




Solution of the equation above will give parameter vector <img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathit{\mathbf{b}}"/>.



<img src="https://latex.codecogs.com/gif.latex?\left\lbrack&space;\begin{array}{cccc}&space;a_f&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;a_f&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;a_f&space;\left\lbrack&space;-\left(q-1\right)\right\rbrack&space;\\&space;a_f&space;\left\lbrack&space;1\right\rbrack&space;&space;&&space;a_f&space;\left\lbrack&space;0\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;a_f&space;\left\lbrack&space;-\left(q-2\right)\right\rbrack&space;\\&space;\vdots&space;&space;&&space;\vdots&space;&space;&&space;\ddots&space;&space;&&space;\vdots&space;\\&space;a_f&space;\left\lbrack&space;q-1\right\rbrack&space;&space;&&space;a_f&space;\left\lbrack&space;q-2\right\rbrack&space;&space;&&space;\cdots&space;&space;&&space;a_f&space;\left\lbrack&space;0\right\rbrack&space;&space;\end{array}\right\rbrack&space;\left\lbrack&space;\begin{array}{c}&space;b\left\lbrack&space;1\right\rbrack&space;\\&space;b\left\lbrack&space;2\right\rbrack&space;\\&space;\vdots&space;\\&space;b\left\lbrack&space;q\right\rbrack&space;&space;\end{array}\right\rbrack&space;=-\left\lbrack&space;\begin{array}{c}&space;a_f&space;\left\lbrack&space;1\right\rbrack&space;\\&space;a_f&space;\left\lbrack&space;2\right\rbrack&space;\\&space;\vdots&space;\\&space;a_f&space;\left\lbrack&space;q\right\rbrack&space;&space;\end{array}\right\rbrack"/>

  

<img src="https://latex.codecogs.com/gif.latex?P_X^{\textrm{ARMA}}&space;\left(f\right)=\frac{\left|{{\mathit{\mathbf{e}}}_{\textrm{vect}}^{\mathit{\mathbf{b}}}&space;}^T&space;\mathit{\mathbf{b}}{\left|\right.}^2&space;\right.}{\;\left|{{\mathit{\mathbf{e}}}_{\textrm{vect}}^{\mathit{\mathbf{a}}}&space;}^T&space;\mathit{\mathbf{a}}{\left|\right.}^2&space;\right.}"/>


<img src="https://latex.codecogs.com/gif.latex?{\mathit{\mathbf{e}}}_{\textrm{vect}}^{\mathit{\mathbf{a}}}&space;={\left\lbrack&space;\begin{array}{cccc}&space;1&space;&&space;e^{-\textrm{j2}\pi&space;f\;}&space;&space;&&space;\cdots&space;&space;&&space;e^{-\textrm{j2}\pi&space;\textrm{fp}\;}&space;&space;\end{array}\right\rbrack&space;}^T"/> 


<img src="https://latex.codecogs.com/gif.latex?{\mathit{\mathbf{e}}}_{\textrm{vect}}^b&space;={\left\lbrack&space;\begin{array}{cccc}&space;1&space;&&space;e^{-\textrm{j2}\pi&space;f\;}&space;&space;&&space;\cdots&space;&space;&&space;e^{-\textrm{j2}\pi&space;\textrm{fq}\;}&space;&space;\end{array}\right\rbrack&space;}^T"/>


```matlab:Code
[y_arma,f]=spestm_obj.psd_arma();
y_arma_st=spestm_obj_st.psd_arma();
```


```matlab:Code
semilogy(f,y_arma);
ylabel('PSD(W/Hz^2) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_13.png](readme_images/figure_13.png)


```matlab:Code
image_y=log10(y_arma_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_14.png](readme_images/figure_14.png)

## MUSIC PSD

```matlab:Code
[y_music,f] = spestm_obj.psd_music();
y_music_st=spestm_obj_st.psd_music();
```


```matlab:Code
semilogy(f,y_music);
ylabel('PSD(W/Hz) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_15.png](readme_images/figure_15.png)


```matlab:Code
image_y=log10(y_music_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_16.png](readme_images/figure_16.png)

## Min-Norm PSD

```matlab:Code
[y_minnorm,f]=spestm_obj.psd_minnorm();
y_minnorm_st=spestm_obj_st.psd_minnorm();
```


```matlab:Code
semilogy(f,y_minnorm);
ylabel('PSD(W/Hz) ');
xlabel('Frequency(Hz)');
legend(x_channel_labels)
grid minor;
```


![figure_17.png](readme_images/figure_17.png)


```matlab:Code
image_y=log10(y_minnorm_st);
imagesc(x_stax,y_stax,image_y)
caxis([min(min(image_y)) max(max(image_y))]);
xlabel('Time(sec)')
ylabel('Freq(Hz)')
set(gca,'YDir','normal')
colorbar;
```


![figure_18.png](readme_images/figure_18.png)

# Comparison of All PSDs

```matlab:Code
lw=1.2;
sep='--';

semilogy(f,y_per,sep,'linewidth',lw); hold on;
plot(f,y_bt,sep,'linewidth',lw)
plot(f,y_cap,sep,'linewidth',lw)
plot(f,y_aryw,sep,'linewidth',lw)
plot(f,y_armc,sep,'linewidth',lw)
plot(f,y_arburg,sep,'linewidth',lw)
plot(f,y_arma,sep,'linewidth',lw)
plot(f,y_music,sep,'linewidth',lw)
plot(f,y_minnorm,sep,'linewidth',lw)
ylabel('PSD(W/Hz) ');
xlabel('Frequency(Hz)');
legend('Periodogram', 'Blackman-Tukey', 'Capon', 'AR Yule-Walker', 'AR Modified Covariance' ...
    , 'AR Burg', 'ARMA', 'MUSIC', 'Min-Norm','Location','best')
grid on;
```


![figure_19.png](readme_images/figure_19.png)

