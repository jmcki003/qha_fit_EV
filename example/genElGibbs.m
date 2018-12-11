
clear

temp=0:25:300;
press = 0;

y=thermalExpansion('glycine_alpha.freq','e-el-b86bpbeXDM.dat',temp,press,'results')


function y = thermalExpansion(freq_file,el_energy_file,Temp,Press,dirName)

%Script to compute Free Energy at given Temp and Volumes
%Self-generates Grueneisen parameters and Fvib calcs
%clear

feedBack = '';

%%%%
% First Check that the mandatory parameters were 
% initialized and set some defaults
%%%%
if(isempty(freq_file))
  feedBack = char(feedBack,'Warning! Frequency file was not set! Exiting...');
  y=feedBack;
  return;
end

if(isempty(el_energy_file))
  feedBack = char(feedBack,'Warning! Energy file was not set! Exiting...');
  y=feedBack;
  return;
end

if(isempty(Temp))
  feedBack = char(feedBack,'Warning! Temperature was not set! Assuming 0');
  Temp = 0;
end

if(isempty(Press))
  feedBack = char(feedBack,'Warning! Pressure was not set! Assuming 0');
  Press = 0;
end

if(isempty(dirName))
  feedBack = char(feedBack,'Warning! Directory was not set! Assuming current directory');
  dirName = pwd;
end

currentDir=pwd;
if ~strcmp(currentDir,dirName) & ~exist(dirName,'dir')
  mkdir(dirName);
end
phononDir=sprintf('%s/%s', dirName, 'phonons');
eosDir=sprintf('%s/%s', dirName, 'eosFits');
fvibDir=sprintf('%s/%s', dirName, 'fvib');
freeEnergyDir=sprintf('%s/%s', dirName, 'freeEnergy');
summaryDir=sprintf('%s/%s', dirName, 'summary');

if ~exist(phononDir, 'dir') 
    mkdir(phononDir);
end
if ~exist(eosDir, 'dir') 
    mkdir(eosDir);
end
if ~exist(fvibDir, 'dir') 
    mkdir(fvibDir);
end
if ~exist(freeEnergyDir, 'dir') 
    mkdir(freeEnergyDir);
end
if ~exist(summaryDir, 'dir') 
    mkdir(summaryDir);
end

%%%%
% Now set some variables that will be used throughout the function
%%%%

%test = fopen('resorcinolAlpha.freq');
test = fopen(freq_file);
ref_found = false;
plus_found = false;
minus_found = false;
ref_freqs = {};
plus_freqs = {};
minus_freqs = {};
x_old = 0;
eosSet = 0;
eosType = '';
acceptedEOS = '';
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',1000);

%constants needed for the calculation
h= 6.62607363e-34; %Js
Na= 6.0221367e23; %1/mol
kb= 1.3806488e-23; %J/K
R= kb*Na;  %J/(mol K)
c= 2.99792458e8; %m/s
count_structures = 0;

%
%Step 0: Extract frequency data from file
%

tline = fgetl(test);
while ~feof(test)
    %disp(tline);
    if contains(tline,'Frequencies')
        blah = strsplit(tline);
        ref_vol = str2double(blah(2));
        ref_found = true;
        tline = fgetl(test);
        x_old = 0;
        count_structures = count_structures + 1;
    elseif contains(tline,'PlusVolume')
        blah = strsplit(tline);
        plus_vol = str2double(blah(2));
        plus_found = true;
        tline = fgetl(test);
        x_old = 0;
        count_structures = count_structures + 1;
    elseif contains(tline,'MinusVolume')
        blah = strsplit(tline);
        minus_vol = str2double(blah(2));
        minus_found = true;
        tline = fgetl(test);
        x_old = 0;
        count_structures = count_structures + 1;
    end
    
    x = str2double(tline);
    
    if (ref_found) && (plus_found) && (minus_found)
        minus_freqs = [minus_freqs x];
    elseif (ref_found) && (plus_found) && (~minus_found)
        plus_freqs = [plus_freqs x];
    elseif (ref_found) && (~plus_found) && (~minus_found)
        ref_freqs = [ref_freqs x];        
    end
    
    if(x-x_old > 2000) && (x_old~=0)
        %fprintf('x: %f\tx_old: %f\n',x,x_old);
        count_structures = count_structures + 1;
    end
    
    if(x ~= 0)
        x_old = x;
    end
    
    tline = fgetl(test);
end

count_structures = count_structures / 3; 
%count_structures = 1;
fprintf('count_structures: %i\n',count_structures);
tmpstr = sprintf('Number of Bruellion zones sampled: %i',count_structures);
feedBack = char(feedBack,tmpstr);

%convert from a cell array to a regular array
ref_freqs = cell2mat(ref_freqs);
plus_freqs = cell2mat(plus_freqs);
minus_freqs = cell2mat(minus_freqs);

%get rid of NaN in the array
ref_freqs(isnan(ref_freqs)) = [];
plus_freqs(isnan(plus_freqs)) = [];
minus_freqs(isnan(minus_freqs)) = [];
fclose(test);

%Now ensure there are no imaginary frequencies
%If so, zero those and the cooresponding modes in 
%the reference frequency lists
ref_index=find(ref_freqs(1:end)<0);
plus_index=find(plus_freqs(1:end)<0);
minus_index=find(minus_freqs(1:end)<0);

if(~isempty(ref_index));
  feedBack = char(feedBack,'Warning! Found an imaginary reference frequency. Zeroing...');
  ref_freqs(ref_index) = 0;
  plus_freqs(ref_index) = 0;
  minus_freqs(ref_index) = 0;
end
if(~isempty(plus_index));
  feedBack = char(feedBack,'Warning! Found an imaginary plus frequency. Zeroing...');
  ref_freqs(plus_index) = 0;
  plus_freqs(plus_index) = 0;
  minus_freqs(plus_index) = 0;
end 
if(~isempty(minus_index));
  feedBack = char(feedBack,'Warning! Found an imaginary minus frequency. Zeroing...');
  ref_freqs(minus_index) = 0;
  plus_freqs(minus_index) = 0;
  minus_freqs(minus_index) = 0;
end


%
%Step 1: Extract electronic energy data from file.
%

%e_el = importdata('e-el-dft-alpha.dat');
e_el = importdata(el_energy_file);
%Volumes = e_el(:,1)';
%Energies = e_el(:,2)';
Volumes = e_el.data(:,1)';
Energies = e_el.data(:,2)';
[minEnergy,pos] = min(Energies);
minVolume = Volumes(pos);

%
%Step 2: Perform EOS fit to the El-Energy surface
%Note that this will not be a true EOS fit hence
%The bulk modulus and its dimensionless derivative
%will not be reported
%

%Murnaghan
F_M = @(a,x) a(1) + (a(2) * x)/a(3) .* (((a(4)./x).^(a(3)))/(a(3) - 1) + 1 ) - ((a(2) * a(4) )/(a(3) - 1));

%Birch-Murnaghan
F_BM = @(a,x) a(1) + (9 * a(2) * a(4))/16 * (((((a(4)./x).^(2/3)-1).^3)*a(3)) + ((((a(4)./x).^(2/3))-1).^2).*(6-4*(a(4)./x).^(2/3)));

%Poirier-Tarantola logarithmic
F_PT = @(a,x) a(1) + (a(2)*a(4)/2)*(log(a(4)./x)).^2 + (a(2)*a(4)/6)*((log(a(4)./x)).^3)*(a(3)-2);

%Vinet
F_V = @(a,x) a(1) + (2*a(2)*a(4)/(a(3)-1)^2)*(2-(5+3*((x/a(4)).^(1/3))*(a(3)-1) - 3*a(3))).*exp(-(3/2)*(a(3)-1).*(((x/a(4)).^(1/3))-1));
    
% Set initial guess parameters and perform least squares fit on the E(V)
% data.
a0 = [minEnergy, 5.0, 8.0, minVolume];
lb = [minEnergy-100000, 1,1, minVolume-50];
ub = [minEnergy+100000, 100,100, minVolume+100];

%fprintf('Trying Murnaghan Fit Parameters:\n');
%Energy_fit_M = lsqcurvefit(F_M,a0,Volumes,Energies,lb,ub,options);
%fprintf('Trying Birch-Murnaghan Fit Parameters:\n');
%Energy_fit_BM = lsqcurvefit(F_BM,a0,Volumes,Energies,lb,ub,options);
%fprintf('Trying Poirier-Tarantola Fit Parameters:\n');
%Energy_fit_PT = lsqcurvefit(F_PT,a0,Volumes,Energies,lb,ub,options);
%fprintf('Trying Vinet Fit Parameters:\n');
%Energy_fit_V = lsqcurvefit(F_V,a0,Volumes,Energies,lb,ub,options);

%fprintf('Trying Murnaghan Fit Parameters:\n');
Energy_fit_M = lsqcurvefit(F_M,a0,Volumes,Energies);
%fprintf('Trying Birch-Murnaghan Fit Parameters:\n');
Energy_fit_BM = lsqcurvefit(F_BM,a0,Volumes,Energies);
%fprintf('Trying Poirier-Tarantola Fit Parameters:\n');
Energy_fit_PT = lsqcurvefit(F_PT,a0,Volumes,Energies);
%fprintf('Trying Vinet Fit Parameters:\n');
Energy_fit_V = lsqcurvefit(F_V,a0,Volumes,Energies);

%fprintf('Murnaghan Fit Parameters:\n');
%fprintf('   A0 = %.5f\n',Energy_fit_M(1));
%fprintf('   V0 = %.3f\n',Energy_fit_M(4));
%fprintf('   B0 = %.3f\n',Energy_fit_M(2));
%fprintf('   B0prime = %.3f\n\n',Energy_fit_M(3));
%Emin_M = Energy_fit_M(1);

%fprintf('Birch-Murnaghan Fit Parameters:\n');
%fprintf('   A0 = %.5f\n',Energy_fit_BM(1));
%fprintf('   V0 = %.3f\n',Energy_fit_BM(4));
%fprintf('   B0 = %.3f\n',Energy_fit_BM(2));
%fprintf('   B0prime = %.3f\n\n',Energy_fit_BM(3));
%Emin_BM = Energy_fit_BM(1);

%fprintf('Poirier-Tarantola Fit Parameters:\n');
%fprintf('   A0 = %.5f\n',Energy_fit_PT(1));
%fprintf('   V0 = %.3f\n',Energy_fit_PT(4));
%fprintf('   B0 = %.3f\n',Energy_fit_PT(2));
%fprintf('   B0prime = %.3f\n\n',Energy_fit_PT(3));
%Emin_PT = Energy_fit_PT(1);

%fprintf('Vinet E(V) Fit Parameters:\n');
%fprintf('   A0 = %.5f\n',Energy_fit_V(1));
%fprintf('   V0 = %.3f\n',Energy_fit_V(4));
%fprintf('   B0 = %.3f\n',Energy_fit_V(2));
%fprintf('   B0prime = %.3f\n\n',Energy_fit_V(3));
%Emin_V = Energy_fit_V(1);
    
Emin_M = Energy_fit_M(1);
Emin_BM = Energy_fit_BM(1);
Emin_PT = Energy_fit_PT(1);
Emin_V = Energy_fit_V(1);
    
nameEOS = sprintf( 'eosFit.txt');
fitEq = fopen(nameEOS,'w');
fprintf(fitEq,'Murnaghan Fit Parameters:\n');
fprintf(fitEq,'   A0 = %.5f\n',Energy_fit_M(1));
fprintf(fitEq,'   V0 = %.3f\n',Energy_fit_M(4));
fprintf(fitEq,'   B0 = %.3f\n',Energy_fit_M(2));
fprintf(fitEq,'   B0prime = %.3f\n\n',Energy_fit_M(3));

fprintf(fitEq,'Birch-Murnaghan Fit Parameters:\n');
fprintf(fitEq,'   A0 = %.5f\n',Energy_fit_BM(1));
fprintf(fitEq,'   V0 = %.3f\n',Energy_fit_BM(4));
fprintf(fitEq,'   B0 = %.3f\n',Energy_fit_BM(2));
fprintf(fitEq,'   B0prime = %.3f\n\n',Energy_fit_BM(3));

fprintf(fitEq,'Poirier-Tarantola Fit Parameters:\n');
fprintf(fitEq,'   A0 = %.5f\n',Energy_fit_PT(1));
fprintf(fitEq,'   V0 = %.3f\n',Energy_fit_PT(4));
fprintf(fitEq,'   B0 = %.3f\n',Energy_fit_PT(2));
fprintf(fitEq,'   B0prime = %.3f\n\n',Energy_fit_PT(3));

fprintf(fitEq,'Vinet E(V) Fit Parameters:\n');
fprintf(fitEq,'   A0 = %.5f\n',Energy_fit_V(1));
fprintf(fitEq,'   V0 = %.3f\n',Energy_fit_V(4));
fprintf(fitEq,'   B0 = %.3f\n',Energy_fit_V(2));
fprintf(fitEq,'   B0prime = %.3f\n\n',Energy_fit_V(3));

fclose(fitEq);
movefile(nameEOS,eosDir);
    
%Now to get some basic statistics
format long
    
errorM=Energies-F_M(Energy_fit_M,Volumes);
errorBM=Energies-F_BM(Energy_fit_BM,Volumes);
errorPT=Energies-F_PT(Energy_fit_PT,Volumes);
errorV=Energies-F_V(Energy_fit_V,Volumes);

avg_error_M=mean(abs(errorM));
avg_error_BM=mean(abs(errorBM));
avg_error_PT=mean(abs(errorPT));
avg_error_V=mean(abs(errorV));

max_error_M=max(abs(errorM));
max_error_BM=max(abs(errorBM));
max_error_PT=max(abs(errorPT));
max_error_V=max(abs(errorV));

name1 = sprintf('murnaghanEOSfitErrors.txt');
mfitEq = fopen(name1,'w');
name2 = sprintf('birch-murnaghanEOSfitErrors.txt');
bmfitEq = fopen(name2,'w');
name3 = sprintf('poirier-tarantolaEOSfitErrors.txt');
ptfitEq = fopen(name3,'w');
name4 = sprintf('vinetEOSfitErrors.txt');
vfitEq = fopen(name4,'w');
name6 = sprintf('allFitErrors.txt');
allfitEq = fopen(name6,'w');
name7 = sprintf('allRelativeFitErrors.txt');
allRelfitEq = fopen(name7,'w');

fprintf(mfitEq,'#Volume(Ang^3) E(kJ/mol) E_fit(kJ/mol) \n');
fprintf(bmfitEq,'#Volume(Ang^3) E(kJ/mol) E_fit(kJ/mol) \n');
fprintf(ptfitEq,'#Volume(Ang^3) E(kJ/mol) E_fit(kJ/mol) \n');
fprintf(vfitEq,'#Volume(Ang^3) E(kJ/mol) E_fit(kJ/mol) \n');
fprintf(allfitEq,'#                True Energy   Murnaghan        Birch-Murnaghan   Poirier-Tarantola\n');%    Vinet           Spline\n');
fprintf(allfitEq,'#Volume(Ang^3)   E(kJ/mol)     E_fit(kJ/mol)    E_fit(kJ/mol)     E_fit(kJ/mol)    \n');%    E_fit(kJ/mol)   E_fit(kJ/mol)\n');
fprintf(allRelfitEq,'#                True Energy   \tMurnaghan        Birch-Murnaghan   Poirier-Tarantola  \n');%  Vinet           Spline\n');
fprintf(allRelfitEq,'#Volume(Ang^3)   E(kJ/mol)     \tE_fit(kJ/mol)    E_fit(kJ/mol)     E_fit(kJ/mol)  \n');%      E_fit(kJ/mol)   E_fit(kJ/mol)\n');
for i=1:size(Volumes')
  fprintf(mfitEq,'%f %f %f\n',Volumes(i),Energies(i),F_M(Energy_fit_M,Volumes(i)));
  fprintf(bmfitEq,'%f %f %f\n',Volumes(i),Energies(i),F_BM(Energy_fit_BM,Volumes(i)));
  fprintf(ptfitEq,'%f %f %f\n',Volumes(i),Energies(i),F_PT(Energy_fit_PT,Volumes(i)));
  fprintf(vfitEq,'%f %f %f\n',Volumes(i),Energies(i),F_V(Energy_fit_V,Volumes(i)));
  fprintf(allfitEq,'%f\t%f\t\t%f\t%f\t%f\n',Volumes(i),Energies(i),F_M(Energy_fit_M,Volumes(i)),F_BM(Energy_fit_BM,Volumes(i)),F_PT(Energy_fit_PT,Volumes(i)));
  fprintf(allRelfitEq,'%f\t%f\t\t%f\t%f\t%f\n',Volumes(i),Energies(i),F_M(Energy_fit_M,Volumes(i))-Energies(i),F_BM(Energy_fit_BM,Volumes(i))-Energies(i),F_PT(Energy_fit_PT,Volumes(i))-Energies(i));
end
fclose(mfitEq);
fclose(bmfitEq);
fclose(ptfitEq);
fclose(vfitEq);
fclose(allfitEq);
fclose(allRelfitEq);

movefile(name1,eosDir);
movefile(name2,eosDir);
movefile(name3,eosDir);
movefile(name4,eosDir);
movefile(name6,eosDir);
movefile(name7,eosDir);


name5 = sprintf( 'statsEOS.txt');
stats = fopen(name5,'w');

%EOSlist=[abs(avg_error_M) abs(avg_error_BM) abs(avg_error_PT) abs(avg_error_V) abs(avg_error_Spline)];
%fullEOSlist=[abs(avg_error_M) abs(avg_error_BM) abs(avg_error_PT) abs(avg_error_V)  abs(avg_error_Spline); abs(max_error_M) abs(max_error_BM) abs(max_error_PT) abs(max_error_V) abs(max_error_Spline)];
%sortedEOSlist=sort(EOSlist);
    
%Now to determine which EOS fit to use

%EOSlist=[avg_error_M avg_error_BM avg_error_PT avg_error_V];
%fullEOSlist=[avg_error_M avg_error_BM avg_error_PT avg_error_V; max_error_M max_error_BM max_error_PT max_error_V];
EOSlist=[avg_error_M avg_error_BM avg_error_PT];
fullEOSlist=[avg_error_M avg_error_BM avg_error_PT; max_error_M max_error_BM max_error_PT];
sortedEOSlist=sort(EOSlist);    
    
chooseM=0;
chooseBM=0;
choosePT=0;
chooseV=0;
choice1=[];

if((abs(avg_error_M) == sortedEOSlist(1)) | (abs(avg_error_M) == sortedEOSlist(2)) )
  choice1=1;
end
if( (abs(avg_error_BM) == sortedEOSlist(1)) | (abs(avg_error_BM) == sortedEOSlist(2)) )
  if(isempty(choice1))
    choice1=2;
  else
    choice2=2;
  end
end
if( (abs(avg_error_PT) == sortedEOSlist(1)) | (abs(avg_error_PT) == sortedEOSlist(2)) )
  if(isempty(choice1))
    choice1=3;
  else
    choice2=3;
  end
end
%if( (abs(avg_error_V) == sortedEOSlist(1)) | (abs(avg_error_V) == sortedEOSlist(2)) )
%  if(isempty(choice1))
%    choice1=4;
%  else
%    choice2=4;
%  end
%end
    

if(fullEOSlist(2,choice1)<fullEOSlist(2,choice2))
  if(choice1==1)
    chooseM=1;
  elseif(choice1==2)
    chooseBM=1;
  elseif(choice1==3)
    choosePT=1;  
  %else(choice1==4)
  %  chooseV=1;
  end
else
  if(choice2==2)
    chooseBM=1;
  elseif(choice2==3)
    choosePT=1;  
  %else
  %  chooseV=1;
  end
end
    
if(chooseM)
  acceptedEOS='Choosing the Murnaghan EOS fit\n';
  fprintf('Choosing the Murnaghan EOS fit\n');
  Energy_fit=Energy_fit_M;
  Emin=Emin_M;
  F = @(a,x) F_M(a,x);
  eosType = 'Murnaghan';
elseif(chooseBM)
  acceptedEOS='Choosing the Birch-Murnaghan EOS fit\n'; 
  fprintf('Choosing the Birch-Murnaghan EOS fit\n');
  Energy_fit=Energy_fit_BM;
  Emin=Emin_BM;
  F = @(a,x) F_BM(a,x);
  eosType = 'Birch-Murnaghan';
elseif (choosePT)
  acceptedEOS='Choosing the Poirier-Tarantola logarithmic EOS fit\n';
  fprintf('Choosing the Poirier-Tarantola logarithmic EOS fit\n');
  Energy_fit=Energy_fit_PT;
  Emin=Emin_PT;
  F = @(a,x) F_PT(a,x);
  eosType = 'Poirier-Tarantola';
%else
%  acceptedEOS='Choosing the Vinet EOS fit\n';    
%  fprintf('Choosing the Vinet EOS fit\n');
%  Energy_fit=Energy_fit_V;
%  Emin=Emin_V;
%  F = @(a,x) F_V(a,x);
%  eosType = 'Vinet'; 
end

fprintf(stats,'%s\n',acceptedEOS);
fprintf(stats,'#E_fit_type             Max Error (kJ/mol)        Average Error (kJ/mol) \n');
fprintf(stats,'Murnaghan\t\t%f\t%f\n',max_error_M,avg_error_M);
fprintf(stats,'Birch-Murnaghan\t\t%f\t%f\n',max_error_BM,avg_error_BM);
fprintf(stats,'Poirier-Tarantola\t%f\t%f\n',max_error_PT,avg_error_PT);
%fprintf(stats,'Vinet\t\t\t%f\t%f\n',max_error_V,avg_error_V);
fclose(stats);
movefile(name5,eosDir);
    
% First the raw E(V) data points
%scatter(Volumes,Energies-minEnergy); 
%hold on; % keep plot active while we add curves
% Now add the E(V) fit
%plot(Volumes,F_M(Energy_fit_M,Volumes)-Emin_M);
%plot(Volumes,F_BM(Energy_fit_BM,Volumes)-Emin_BM);
%plot(Volumes,F_PT(Energy_fit_PT,Volumes)-Emin_PT);
%plot(Volumes,F_V(Energy_fit_V,Volumes)-Emin_V);
%xlabel('Volume (Ang^3)');
%ylabel('Energy (kJ/mol)');
%legend('Raw Data','M','BM','PT');
%legend('Raw Data','M','BM','PT','V');
%hold off;

%
%Step 3: Generate the Grueneisen Parameters
%

%grun = -1.*((log(plus_freqs) - log(minus_freqs))./(log(plus_vol) - log(minus_vol)));
grun = -1.*((log(plus_freqs./minus_freqs))./(log(plus_vol./minus_vol)));
grun(isnan(grun)) = 0;

name = sprintf('grueneisenParams.txt');
grunGen = fopen(name,'w');
fprintf(grunGen,'Grueneisen Parameters\n');
for i = 1:1:size(grun(:))
    fprintf(grunGen,'%f\n',grun(i));
end
fclose(grunGen);
movefile(name,dirName);

%
%Step 4: Now establish the Volumes, Pressures, and Temperatures we wish to
%calculate over
%

%Volume to extrapolate from
%v_min = Energy_fit(4)-100;
%v_max = Energy_fit(4)+100;
%vol_array = v_min:v_max;
v_min = min(Volumes(:));
v_max = max(Volumes(:));
vol_array=Volumes;

%Pressure in GPa
for press = Press
  holdName = sprintf('final_resultsP%1.2f.txt', press);
  results = fopen(holdName,'w');
  fprintf(results,'#Temperature (K)     Optimal Volume (Ang^3)    Free Energy (kJ/mol)    Entropy (J/(mol K))    Enthalpy (kJ/mol)    Fvib (kJ/mol)    Internal Energy (kJ/mol)\n');
  
  %Loop over a set of temperatures
  for temp = Temp
    %Reset Fvib values to an empty array
    fprintf('\n\nTemp = %f\n',temp);
    fvib = {};
    hvib = {};
    svib = {};
    cv = {};
    
    %count_structures'
    %vol_array'
    %loop over a set of volumes
    for vol = vol_array

      %calculate new frequency
      %units: m/s * (cm/m) * (1/cm) = 1/s = Hz
      w =  ( c * 100 ) .* ref_freqs .* (( vol / ref_vol ).^( -1 .* grun ));

     % vol
     % w'
    %nameTempor = sprintf( 'fvibEnergyT%iP%1.2fV%1.2f.txt', temp, press,vol );
    %fileID = fopen(nameTempor,'w');
    %fprintf(fileID,'#Volume (Ang^3)   Helmholtz Free Energy (kJ/mol)\n');
    %for i = w
    %  fprintf(fileID,'%f\n',w);
    %end
    %fclose(fileID);
    %movefile(nameTempor,phononDir);
    
    
      %Calculate zero point energy contribution
      zpe = (w .* h) ./ 2;
      zpe = sum(zpe);
    
      %If relevant calculate other term
      if temp ~= 0
       %units: J/K * K = J
       h_vib_term2 = (h .* w) ./ ( exp( ( h .* w )./( kb * temp )) - 1);
       h_vib_term2 (~isfinite(h_vib_term2))=0;
       h_vib_term2 = sum(h_vib_term2);
       s_vib = ((h .* w) ./(temp .* ( exp( ( h .* w )./( kb * temp )) - 1))) - (( kb ) .* log( 1 - exp( ( -1 * h .* w )./( kb * temp ))));
       s_vib(~isfinite(s_vib))=0;
       s_vib= sum(s_vib);   
       Cv = (Na/(count_structures*kb)).*( ( ((h .* w) ./(temp .* ( exp( ( h .* w )./( kb * temp )) - 1))).^2).* exp( (h .* w )./( kb * temp )));
       Cv(~isfinite(Cv))=0;
       Cv = sum(Cv);     
       %fprintf('entr_enth: %f kJ/mol\n',entropy_enthalpy);       
      else
        entropy_enthalpy = 0.0;
        s_vib = 0;
        h_vib_term2 = 0;
        Cv = 0;
      end
      
      %Sum to create Fvib
      %units: (1/mol) * (J + J) * (kJ / J) = kJ/mol
      s_vib =  Na * (s_vib) / (1000 * count_structures);
      h_vib =  Na * (zpe + h_vib_term2) / (1000 * count_structures);
      f_vib = h_vib - (temp * s_vib);
      fvib = [fvib f_vib];
      hvib = [hvib h_vib];   
      svib = [svib s_vib];
      cv = [cv Cv];

      %fprintf('Temp: %f kJ/mol  Vol : %f Ang^3\n',temp,vol);   
      %fprintf('fvib: %f kJ/mol\n',f_vib);
    end
    
    %Call to calculate the gibbs free energy
    fvib = cell2mat(fvib);
    hvib = cell2mat(hvib);
    svib = cell2mat(svib);
    cv = cell2mat(cv);
    
    %y = 'pause'
    %return 
    fvib_spline = spline(vol_array,fvib);
    spline_fvib_func = @(v) ppval(fvib_spline,v);
  
    hvib_spline = spline(vol_array,hvib);
    spline_hvib_func = @(v) ppval(hvib_spline,v);  
  
    svib_spline = spline(vol_array,svib);
    spline_svib_func = @(v) ppval(svib_spline,v); 
  
    cv_spline = spline(vol_array,cv);
    spline_cv_func = @(v) ppval(cv_spline,v); 
    
    %pv_spline = spline(vol_array,pv);
    %spline_pv_func = @(v) ppval(pv_spline,v); 
    
    %hvib'
    %Energies'
    %pv'
    %gibbs'
    %gibbs = fvib + Energies + pv;
    
    %Free_energy = @(v) F(Energy_fit,v) + f_vib_func(v);
    PVterm = @(v) press * v * (1.0e-24) * Na; 
    Free_energy = @(v) F(Energy_fit,v) + spline_fvib_func(v) + PVterm(v); 
    
    %Free_energy = @(v) F(Energy_fit,v);
    
    %Now search over all available volumes for the optimal Gibbs' value
    [Vmin, Gmin, info, output] = fminbnd(Free_energy, v_min, v_max);
    
    % Print out the optimal volume and free energy
    fprintf('Optimal Volume and Free Energy:\n%f Ang^3\t%f kJ/mol\n',Vmin,Gmin);
    Smin = spline_svib_func(Vmin);
    Hmin = spline_hvib_func(Vmin);
    Fvibmin = spline_fvib_func(Vmin);
    CVmin = spline_cv_func(Vmin);
    PVmin = press * Vmin * (1.0e-24) * Na;
    %EelMin = Gmin - PVmin - Fvibmin;
    EelMin = F(Energy_fit,Vmin);
    fprintf('Optimal Entropy, Enthalpy and Cv:\n%f J/mol\t%f kJ/mol\t%f J/(mol K)\n',Smin*1000,Hmin,CVmin);
    fprintf('PV term :\t%f kJ/mol,\tE_el = \t%f\n',PVmin,EelMin);
  
    %
    %Step 5: Store results in 'logical' files
    %
  
    name = sprintf( 'freeEnergyT%iP%1.2f.txt', temp, press );
    fileID = fopen(name,'w');
    fprintf(fileID,'#Volume (Ang^3)   Free Energy (kJ/mol)\n');
    for i = vol_array
      fprintf(fileID,'%f %f\n',i,Free_energy(i));
      %fprintf('%f %f\n',i,Free_energy(i));
    end

    fclose(fileID);
    movefile(name,freeEnergyDir);
    
    name = sprintf( 'fvibEnergyT%iP%1.2f.txt', temp, press );
    fileID = fopen(name,'w');
    fprintf(fileID,'#Volume (Ang^3)   Helmholtz Free Energy (kJ/mol)\n');
    for i = vol_array
      fprintf(fileID,'%f %f\n',i,spline_fvib_func(i));
    end
    fclose(fileID);
    movefile(name,fvibDir);
    
    %fprintf(results,"%f   %f   %f   %f   %f   %f     %f\n",temp,Vmin,Fmin,Smin*1000,Hmin,CVmin,F(Energy_fit,Vmin));
    fprintf(results,"%f   %f   %f   %f   %f   %f     %f\n",temp,Vmin,Gmin,Smin*1000,Hmin,Fvibmin,EelMin);
    
    wMin =  ( c * 100 ) .* ref_freqs .* (( Vmin / ref_vol ).^( -1 .* grun ));

    namePhonon = sprintf( 'phononsT%iP%1.2fV%1.2f.txt', temp, press, Vmin );
    phonon = fopen(namePhonon,'w');
    fprintf(phonon,'Phonons at Temperature %i K Pressure %1.2f GPa and Volume  %1.2f Ang^3\n', temp, press, Vmin );
    for i = 1:1:size(wMin(:))
      fprintf(phonon,'%f\n',wMin(i));
    end
    fclose(phonon);
    movefile(namePhonon,phononDir);

  end
  
  fclose(results);
  movefile(holdName,summaryDir);
end

disp('Made it to the end')
feedBack = char(feedBack,'Job Complete!');
y = feedBack;
return
end
