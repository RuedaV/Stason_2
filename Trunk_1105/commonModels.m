filename = {'CAC40.csv', 'FTSE100.csv', 'HandSeng.csv', ... 
    'Nikkei225.csv', 'S&P500.csv'};    
datasource = xlsread(filename{3});
n = length(datasource(:,1));
d = datasource(1:n, 1);
r = datasource(1:n, 2);
disp (filename{fileIndex});
model = gjr('GARCHLags',1,'ARCHLags',1,'LeverageLags',1);
EstMdl = estimate(model,r);