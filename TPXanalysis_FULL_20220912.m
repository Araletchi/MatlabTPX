clc; clear;

% **********************************************************************************************************************************
%
% DATA STRUCTURE INFO
%
%   {T}: events splitted 1 by 1 including  x_pos  #  y_pos  #  ToA  #  ToT  #  Blob_label  #  [Blank]  #  Corrected_x  #  Corrected_y  #  Inner_pixels 
%   {data}: pixel info from each TPX3 file (not splitted)
%   {IM}: Matrix representation of each TPX3 file
%   {labels}: IM matrices labelling each event
%   {M}: same as {data} without event labels
%   {props}: properties of each TPX3 file regions of interest (pixels, size, centroids, ...)
%   {tracks}: events resulting from each TPX3 file (same as {T} but for each non-empty TPX3 file)
%
%   **** Pixel size calculation only applies for cricket noise measurements
%
%
% OUTPUT: 
%   date_time-Resulting_tracks.txt file shows track information using CERF algorith and track shape
%   date_time-TIFF_files.txt shows TIFF file info (number of pixels and intensity)
%   date_time_Resulting_image.txt shows the merged image of all tracks after TIFF files
%
% **********************************************************************************************************************************

[file_sel,path1]=uigetfile('*.*', 'Select file: ');                      % Select one file from folder
[filepath,name,ext]=fileparts(file_sel);  
w=warning('off','all');

filelist=fullfile(path1,append('*',ext));                                   % Filelist with the selected file extension
source=dir(filelist);

tpfiles=dir([path1 '\*.tpx3']);
tifiles=dir([path1 '\*.tiff']);

ltpx=length(tpfiles);
ltif=length(tifiles);


% ******************** (1) TPX3 FILES ANALYSIS (IF EXIST) ********************

if ltpx>0

    m={};                                                                   % Cell to store 
    n={};
    omit={};
        
    % Read files, save data as unique file and delete temporary files and variables
    
    for k=1:ltpx                                                            % Reads all the files one by one
        imglist=tpfiles(k).name;
        tpxfile=fullfile(tpfiles(k).folder,imglist);
        OutputFileName=append(datestr(now,'yyyymmdd-'),tpfiles(k).name,'.txt');  % Output TXT to be named after input file # 'yyyymmdd_HHMM-' optional
    
        tpx2txt(tpxfile,OutputFileName);                                    % Run the function
        data=importdata(OutputFileName);                                    % Extracts values from output file      
        
        files=dir('2022*.txt');                                             % Remove empty TXT files
        for k=1:numel(files)
            if files(k).bytes==0
                delete(files(k).name);
            end
        end
        
        files2=dir('2022*.txt'); 
        fileout=append(datestr(now,'mmdd-'),'TPXcricket','.txt');
        fout=fopen(fileout,'wt');
        
        for j=1:length(files2)
            fin=fopen(files2(j).name);
            temp=fread(fin,'uint8');
            fprintf(fout,'%s\n',files2(j).name,'',temp);
            fclose(fin);
        end
    
        if isempty(data)== 0                                                % If TXT file is not empty, reads content
            xpix=data(:,1);                                                 % Stores variables for plot
            ypix=data(:,2);
            toa=data(:,3);
            tot=data(:,4);
              
            n{k}=[xpix ypix toa tot];                                       % Cell to store the x,y values of each TPX3 file             
        end                      
    end
    
    m=n(~cellfun('isempty',n));                                             % Remove the empty cells (if exist)
    
    files=dir('2022*.txt');                                                 % Remove the TXT files after merging in a single file
        for k=1:numel(files)
            delete(files(k).name);
        end
end    
    
    
% ************************** (2) CREATE CELL TO SPLIT EACH FILE VALUES **************************
    
lem=length(m);
    
% Omit one-line values (single pixels) for the analysis and assign the ToT value to the (x,y) corresponding pairs
    
    for k=1:lem                                                                  
        meas=size(m{k});
        
        if meas(1)>1                                                        % Omit 1-line values
            omit{k}=m{k};
        else
            omit{k}=[];
        end
    end

M=omit(~cellfun('isempty',omit));                                           % Remove the empty cells (if exist)
    
clearvars ans data ext file_sel filelist fileout filepath files files2 fin fout imglist j k lem ltpx m meas n name omit 
clearvars OutputFileName path1 source temp toa tot tpfiles tpxfile w xpix ypix  


% ******************** (3) IDENTIFY CLUSTERS ******************** 

data={};

for i=1:length(M)
    data{1,i}=[(M{1,i}(:,1)+1) (M{1,i}(:,2)+1) M{1,i}(:,3) M{1,i}(:,4)];    % TPX (x,y) values increased in 1 to fit matrix indices
end

b={};
IM={};
labels={};

for i=1:length(data)
    valuesx=data{1,i}(:,1);
    valuesy=data{1,i}(:,2);

    b{1,i}=zeros(256);
        
    for j=1:length(data{1,i})
        b{1,i}(valuesx(j),valuesy(j))=256;                                  % Fill each 256x256 matrix with non-zero values if activated pixel        
    end

    IM{1,i}=uint8(b{1,i});

    [labels{1,i}]=bwlabel(IM{1,i},8);
    props{1,i}=regionprops(labels{1,i},IM{1,i},'all');
end

xval={}; yval={};                                                           % Empty cells to store x,y blobs length


% ************************** (4) SPLIT EVENTS **************************
ld=length(data);

for i=1:ld
    for j=1:length(data{1,i})
        data{1,i}(j,5)=labels{1,i}(data{1,i}(j,1),data{1,i}(j,2));          % Add props label info to data
    end
end


% ************************** (5) SPLIT PAIRS UPON CONDITION (DISTANCE) **************************
tracks={};

for i=1:ld
    data{1,i}=sortrows(data{1,i},[5]);
end


for i=1:length(data)
    if (max(data{1,i}(:,5),[],'all'))==1
        tracks{1,i}=data{1,i};

    else
        a=data{1,i}(:,5);
        for j=1:(max(a))
            if j==1
                tracks{j,i}=data{1,i}(1:(max(find(a==1))),:);
            else
                for j=2:(max(a))              
                    tracks{j,i}=data{1,i}((max(find(a==(j-1)))+1):(max(find(a==(j)))),:); 
                end
            end
        end
    end
end

T=tracks(~cellfun('isempty',tracks));                                       % Remove the empty cells (if exist)                

clearvars a b i j ld valuesx valuesy


% ******************** (6) CLUSTERS IDENTIFICATION ********************

% Create a file to save results
fid=fopen(append(datestr(now,'yyyymmdd_HHMM-'),'Resulting_tracks','.txt'),'wt');

lT=length(T);
Tmeas={};

for k=1:lT
    if length(T{k}(:,1))>1                                                  % Skip single pixels from analysis
        z=[(T{k}(:,1)-min(T{k}(:,1))+1) (T{k}(:,2)-min(T{k}(:,2))+1)];
        A=[];
        
        for i=1:size(z,1)                                                   % If x,y pair value exists, writes 1 (otherwise, is 0)     
            A(z(i,1),z(i,2))=1;
        end
        
        % Save new coordinate values together with the tracks information

        % Transform A to "plot" the track shape
        B=A;                                                                % Create a copy       
        B=num2str(B); B(B=='1')='X';  B(B=='0')='·';                        % Convert to string and replace the 1's with X's
        str=char(B);
        
        dx=size(A,1);   dy=size(A,2);                                       % Matrix size                     
    
        inn=[];
    
        for i=1:length(z)
            z1=z(i,1); z2=z(i,2); 
            
            for j=1:length(z)
                nb1=0;nb2=0;nb3=0;nb4=0;nb5=0;nb6=0;nb7=0;nb8=0;
                if ismember([z1-1 z2-1],z,'rows'); nb1=1; end;
                if ismember([z1-1 z2],z,'rows'); nb2=1; end;
                if ismember([z1-1 z2+1],z,'rows'); nb3=1; end;
                if ismember([z1 z2-1],z,'rows'); nb4=1; end;
                if ismember([z1 z2+1],z,'rows'); nb5=1; end;
                if ismember([z1+1 z2-1],z,'rows'); nb6=1; end;
                if ismember([z1+1 z2],z,'rows'); nb7=1; end;
                if ismember([z1+1 z2+1],z,'rows'); nb8=1; end;
            end
            
            inn(i)=nb1+nb2+nb3+nb4+nb5+nb6+nb7+nb8;
            
            % Save new coordinates x,y values and inner pixels for each track T
            
            T{k}(i,7)=z1;
            T{k}(i,8)=z2;
            T{k}(i,9)=inn(i);

        end
        
        inner=sum(inn>4);                                                   % Find number of pixels with over than 4 neighbors                                                 
    
        n=sum(A,'all');                                                     % Total number of pixels of the track
        
        % Maximum height for the track elements
        for i=1:dy
            w(i)=max(sum(A(:,i),'all'));                
        end
        
        % Track length, width and angle
        L=sqrt(dx^2+dy^2);                                                          
        theta=atand(dy/dx);                                                         
        W=(max(w))*cosd(theta);                                                                 
        
        % Ratio and density values for track identification
        mindim=min(dx,dy);
        maxdim=max(dx,dy);
    
        if (mindim/maxdim<=1) && (mindim/maxdim>=0.8) && (n/(mindim*maxdim)>0.75)
            ratio=mindim/maxdim;
            dens=n/(dx*dy);
        else
            ratio=L/W;                                                          % Ratio 
            dens=n/(L*W);                                                       % Track density 
    
        end
    
        maxis=min(dx,dy);       Maxis=max(dx,dy);                               % Minor (and major) axis size in pixels        
        
            % CERF CLUSTERING ALGORITHM ***********************************************
            % Small blob: n=0 && n=<4
            % Heavy track: inner>4 && ratio>1.25 && dens>0.3
            % Heavy blob: inner>4 && ratio<1.25 && dens>0.5
            % Medium blob: 1<=inner<=4 && ratio<1.25 && dens>0.5
            % Straigth track: ratio>8 && maxis<3
            % Light track: Not straight track
            
            if (inner==0 && n<=4)
                fprintf(fid,'Track %d: Small blob. \n Total pixels: %d - Inner pixels: %d - L/W ratio: %0.2f - Density: %0.2f \n\n',k,n, inner, ratio, dens);
                fprintf(fid,' x    y       ToA      ToT\n --------------------------\n');
                fprintf(fid,'%g  %g  %0.9f  %g \n\n', [T{k}(:,1),T{k}(:,2),T{k}(:,3),T{k}(:,4)].');
                fprintf(fid,'  Track shape:\n');
                for j=1:size(B,1)
                    fprintf(fid,'  %s\n',B(j,:));
                end
                fprintf(fid,'\n\n');
                
            elseif (inner>4 && ratio>1.25 && dens>0.3)
                fprintf(fid,'Track %d: Heavy track. \n Total pixels: %d - Inner pixels: %d - L/W ratio: %0.2f - Density: %0.2f \n\n',k,n, inner, ratio, dens);
                fprintf(fid,' x    y       ToA      ToT\n --------------------------\n');
                fprintf(fid,'%g  %g  %0.9f  %g \n\n', [T{k}(:,1),T{k}(:,2),T{k}(:,3),T{k}(:,4)].');
                fprintf(fid,'  Track shape:\n');
                for j=1:size(B,1)
                    fprintf(fid,'  %s\n',B(j,:));
                end
                fprintf(fid,'\n\n');
    
            elseif (inner>4 && ratio<1.25 && dens>0.5)
                fprintf(fid,'Track %d: Heavy blob. \n Total pixels: %d - Inner pixels: %d - L/W ratio: %0.2f - Density: %0.2f \n\n',k,n, inner, ratio, dens);
                fprintf(fid,' x    y       ToA      ToT\n --------------------------\n');
                fprintf(fid,'%g  %g  %0.9f  %g \n\n', [T{k}(:,1),T{k}(:,2),T{k}(:,3),T{k}(:,4)].');
                fprintf(fid,'  Track shape:\n');
                for j=1:size(B,1)
                    fprintf(fid,'  %s\n',B(j,:));
                end
                fprintf(fid,'\n\n');
    
            elseif (inner>=1 && inner<=4 && ratio<1.25 && dens>0.5)
                fprintf(fid,'Track %d: Medium blob. \n Total pixels: %d - Inner pixels: %d - L/W ratio: %0.2f - Density: %0.2f \n\n',k,n, inner, ratio, dens);
                fprintf(fid,' x    y       ToA      ToT\n --------------------------\n');
                fprintf(fid,'%g  %g  %0.9f  %g \n\n', [T{k}(:,1),T{k}(:,2),T{k}(:,3),T{k}(:,4)].');
                fprintf(fid,'  Track shape:\n');
                for j=1:size(B,1)
                    fprintf(fid,'  %s\n',B(j,:));
                end
                fprintf(fid,'\n\n');
    
            elseif (inner==0 && ratio>8 && maxis<3)
                fprintf(fid,'Track %d: Straight track. \n Total pixels: %d - Inner pixels: %d - L/W ratio: %0.2f - Density: %0.2f \n\n',k,n, inner, ratio, dens);
                fprintf(fid,' x    y       ToA      ToT\n --------------------------\n');
                fprintf(fid,'%g  %g  %0.9f  %g \n\n', [T{k}(:,1),T{k}(:,2),T{k}(:,3),T{k}(:,4)].');
                fprintf(fid,'  Track shape:\n');
                for j=1:size(B,1)
                    fprintf(fid,'  %s\n',B(j,:));
                end
                fprintf(fid,'\n\n');
    
            else 
                fprintf(fid,'Track %d: Light track. \n Total pixels: %d - Inner pixels: %d - L/W ratio: %0.2f - Density: %0.2f \n\n',k,n, inner, ratio, dens);
                fprintf(fid,' x    y       ToA      ToT\n --------------------------\n');
                fprintf(fid,'%g  %g  %0.9f  %g \n\n', [T{k}(:,1),T{k}(:,2),T{k}(:,3),T{k}(:,4)].');
                fprintf(fid,'  Track shape:\n');      
                for j=1:size(B,1)
                    fprintf(fid,'  %s\n',B(j,:));
                end
                fprintf(fid,'\n\n');
               
            end
    end
end
    fclose(fid);

clearvars A ans B dens dx dy fid i j k inn inner L lT maxdim maxis Maxis mindim n
clearvars nb1 nb2 nb3 nb4 nb5 nb6 nb7 nb8 question ratio str theta Tmeas 
clearvars w W xval yval z z1 z2


% ******************** (7) CRICKET NOISE SIZE / EVENTS SIZE ******************** 

question=input('Average pixel size (x,y) for noise cricket size calculation? Y/N [Y]:','s');

if question=='Y'

    for i=1:length(props)
        for j=1:length(props{1,i})
            xval{i,j}=length(props{1,i}(j).Image(:,1));
            yval{i,j}=length(props{1,i}(j).Image(1,:));
        end
    end
        x_size=mean([xval{:}]);
        y_size=mean([yval{:}]);
        
        fprintf('Mean x cricket noise size: %0.2f \nMean y cricket noise size: %0.2f \n', x_size, y_size);
        
        clearvars i j xval yval
    else
        fprintf(' ');
end


% ******************** (8)  TIFF FILES ANALYSIS (IF EXIST) ********************

if ltif>0

    R=uint16(zeros(256));                                                       % Initiliaze matrix to save resulting image
    
    % Create a file to store the non-empty images list
    fileID=fopen(append(datestr(now,'mmdd_HHMM-'),'TIFF_files','.txt'),'wt');
    
    for k=1:ltif                                                               % Along all the files
        imglist=tifiles(k).name;
        allimg=fullfile(tifiles(k).folder,imglist);
        A=imread(allimg);
        a=sum(A,'all');
        b=length(find(A(:)>1));                                                 % Activated pixels
    
        if b<=1                                                                 % Avoid single pixels hotspots                                                
            R=R;                    
        else %(a>0 && max_a<15000)                                              % Saves to a single matrix and non-zero values to a TXT file
            R=R+A;
            fprintf(fileID, 'File: %s  - Sum: %d Pixels: %d \n', imglist,a,length(find(A(:)>1)));
        end
    end
    
    fclose(fileID);
    
    if sum(R,'all')~=0
        figure(1)
        result=imagesc(R);                                                          % Show image
        title('Resulting image (TIFF composition)');
        colorbar;
        % Save resulting image including date
        hgexport(gcf,append(datestr(now,'mmdd_HHMM-'),'Resulting_image'), hgexport('factorystyle'), 'Format', regexprep('tiff','[.]',''));
    end
end


% ******************** (99) FUNCTION FOR TPX3 FILES ********************

function tpx2txt(tpx3file, fileout)
  D=dir(tpx3file);
  filesize=D.bytes;                           
  [filepath,name,ext]=fileparts(tpx3file); 

  if filesize==0                                                            % Check file integrity, exit if size=0
    fprintf('Sum Tin Wong \n\n');
    return
  end

  fid=fopen(tpx3file,'rb');              
      
  data=typecast(uint32(fread(fid,filesize,"uint32",0,"ieee-le")),"uint64"); % Read all data from the file into array according to format and close the file
  fclose(fid);

  fprintf('\nFile: %s%s  \nx_pix|y_pix|    ToA (ns)    | ToT \n --------------------------------- \n',name,ext);

  ldata=uint64(length(data)-1);                                             % Length of data for the while loop 
  dout=strings(ldata,1);                                                    % Convert length to string
  jj=1;
  i=uint64(1);                                                              % Convert the starting ldata to proper format
    
  while ldata>i
    temp64=data(i);                                                         % Analyse the data packet header     
    i=i+1;
    val_T=char(bitand(temp64,uint64(0xFF)));
    val_p=char(bitand(bitshift(temp64,-8),uint64(0xFF)));
    val_x=char(bitand(bitshift(temp64,-16),uint64(0xFF)));
    val_3=char(bitand(bitshift(temp64,-24),uint64(0xFF))); 

    if (val_T=="T" && val_p=="P" && val_x=="X" && val_3=="3")
      chipnr = bitand(bitshift(temp64, -32),uint64(0xFF));
      packets = bitshift(temp64,-51)-1;                                     % Split packets

      for ii = i:(i+packets)                                                % Goes through every packet 
        temp64 = data(ii);
        header = bitshift(temp64, -60);

        switch (header)
          case {uint64(0x6)}                                                % Parse TDC data packet 
            tcount = bitand(bitshift(temp64, -44), 0x0000000000000FFF);
            coarsetime = bitand(bitshift(temp64, -9),0x00000007FFFFFFFF);
            coarse_time = double(coarsetime)/320e6;
            timefine = bitand(bitshift(temp64, -5), 0x000000000000000F);
            TDC_timestamp = coarse_time + double(timefine)*260e-12;
            % dout(jj,1)=sprintf("tdc: %15.12f",TDC_timestamp);             % TDC data to string array, comment to omit the output
            % jj=jj+1;

          case {uint64(0xB)}                                                % Parse pixel hit data packet               
            dcol=bitshift(bitand(temp64,0x0FE0000000000000),-52);           % -60+7+1       % Calculating pixel coordinates
            spix=bitshift(bitand(temp64,0x001F800000000000),-45);           % -53+6+2
            pix0=bitshift(bitand(temp64,0x0000700000000000),-44);           % -47+3
            pix_x=int32(dcol+bitshift(pix0,-2));
            pix_y=int32(spix+bitand(pix0,uint64(0x3)));
                                                
            spidrtime=bitand(temp64,0x000000000000FFFF);                    % Calculating timings                       
            ToA=bitand(bitshift(temp64,-30), 0x0000000000003FFF);
            ToT=bitand(bitshift(temp64,-20), 0x00000000000003FF);
            fToA=bitand(bitshift(temp64,-16), 0x000000000000000F);
            cToA=ToA*25;
            spidrtime_ns=spidrtime*25*16384;

            ToA_ns=ToA*25;                                                  % Time in ns
            ToT_ns = ToT*25;
            chip_ToA=double(spidrtime_ns+cToA-(fToA*25)/16)/1e9;            % Time in ns

        switch(chipnr +3)                                                   % Recalculating pixel coordinates for the quad (2x2) chips
          case {3}

          case {0}
            pix_x=pix_x+260;
            pix_y=pix_y;

          case {1}
            pix_x=255-pix_x+260;
            pix_y=255-pix_y+260;

          case {2}
            pix_x=255-pix_x;
            pix_y=255-pix_y+260; 
        end
      
      % Pixel data to string array
      dout(jj,1)=sprintf("%5d %5d %15.12f %6d",pix_x,pix_y,chip_ToA,ToT_ns);
      jj=jj+1;
      
      fprintf("%5d %5d %15.12f %6d\n",pix_x,pix_y,chip_ToA,ToT_ns);
      % case {uint64(0x4)} % #"global time"

        end
     end
    
     i=i+packets;

  end
end

% Output data
jj=jj-1;
writematrix(dout(1:jj,1),fileout);
%fprintf("%5d %5d %15.12f %6d \n",pix_x,pix_y,chip_ToA,ToT_ns);

end
