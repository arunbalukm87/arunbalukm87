% Open NIRS_KIT and create the mat file

%Select the Hbx_conc.mat for performing correlation.
uiopen("*.mat");

%%
% Calculating Correlation matrix

R_HbO=corr(Conc.HbO);
R_HbR=corr(Conc.HbR);
R_HbT=corr(Conc.HbT);
%%
%calculating Partial Correlation


Rp_HbO=partialcorr(Conc.HbO);
Rp_HbR=partialcorr(Conc.HbR);
Rp_HbT=partialcorr(Conc.HbT);
% splitting the correlation matrix into left and right hemipshere
RpHbO_Left=Rp_HbO(1:10,1:10);
RpHbO_right=Rp_HbO(11:20,11:20);
RpHbR_Left=Rp_HbR(1:10,1:10);
RpHbR_right=Rp_HbR(11:20,11:20);
RpHbT_Left=Rp_HbT(1:10,1:10);
RpHbT_right=Rp_HbT(11:20,11:20);

% removing the diagonal
RpHbO_Left=weight_conversion(RpHbO_Left,'autofix');
RpHbO_right=weight_conversion(RpHbO_right,'autofix');
%analysing only positive correlations
RpHbO_left=max(RpHbO_Left,0);
RpHbO_right=max(RpHbO_right,0);
RpHbR_Left=max(RpHbR_Left,0)
RpHbR_right=max(RpHbR_right,0);
RpHbT_Left=max(RpHbT_Left,0)
RpHbT_right=max(RpHbT_right,0);

%fisher transform
Z_RpHbO_left=atanh(RpHbO_left);
Z_RpHbO_right=atanh(RpHbO_right);
Z_RpHbR_left=atanh(RpHbO_Left);
Z_RpHbR_right=atanh(RpHbR_right);
Z_RpHbT_left=atanh(RpHbT_Left);
Z_RpHbT_right=atanh(RpHbT_right);

%%
tot=numel(RpHbO_left);
val=sum(RpHbO_left(:)==0);
SpR=val/tot;

% RpHbO_LeftZ=weight_conversion(RpHbO_Left,'normalize');
%%
%Plotting the heat map 
    % figure;
    % heatmap(R_HbO,"Colormap",jet,"Title",'HbO Corr Matrix');figure;
    % heatmap(R_HbR,"Colormap",jet,"Title",'HbR Corr Matrix');figure;
    % heatmap(R_HbT,"Colormap",jet,"Title",'HbT Corr Matrix');


figure;
heatmap(Rp_HbO,"Colormap",jet,"Title",'HbO pCorr Matrix');figure;
heatmap(Rp_HbR,"Colormap",jet,"Title",'HbR pCorr Matrix');figure;
heatmap(Rp_HbT,"Colormap",jet,"Title",'HbT pCorr Matrix');

%%
% Fisher transform
    % r_hbo=Rp_HbO(:);    % z=.5.*log((1+r_hbo)./(1-r_hbo));
    % Z_hbo=atanh(r_hbo);
Z_HbO=atanh(Rp_HbO);
Z_HbR=atanh(Rp_HbR);
Z_HbT=atanh(Rp_HbT);
Z_HbO = weight_conversion(Z_HbO,'autofix');
Z_HbR = weight_conversion(Z_HbR,'autofix');
Z_HbT = weight_conversion(Z_HbT,'autofix');
% figure;
% heatmap(Z_HbO,"Colormap",jet,"Title",'HbO pCorr Matrix');figure;
% heatmap(Z_HbR,"Colormap",jet,"Title",'HbR pCorr Matrix');figure;
% heatmap(Z_HbT,"Colormap",jet,"Title",'HbT pCorr Matrix');
% 


%normalise between 0 and 1
Z_HbO_norm = weight_conversion(Z_HbO,'normalize');
Z_HbR_norm = weight_conversion(Z_HbR,'normalize');
Z_HbT_norm = weight_conversion(Z_HbT,'normalize');
%%
fc=Z_HbO_norm;
imagesc(fc); colorbar; %Notice negative weights


%% ROInames
N=20;
pre = 'CH';
names = {};
for k = 1:20
    names = [names;strcat([pre,num2str(k,'%02d')])];
end
%% COmmunity structure/Modularity
  
    for i=1:100
                                 [Ci0 Q0]= community_louvain(fc,[],[],'negative_sym');
                                 Ci_iter(i,:)= Ci0;
                          end
                          Ci=transpose(Ci_iter);

                          D= agreement(Ci);
                          
                          ci0=consensus_und(D,80,100);
 %% Visualize Modular Structure of Adjacency Matrix %%
        figure;
        [X,Y,INDSORT] = grid_communities(Ci0); % call function
        imagesc(fc(INDSORT,INDSORT));           % plot ordered adjacency matrix
        hold on;                                 % hold on to overlay community visualization
        
        plot(X,Y,'y','linewidth',3);
        title('Community Structure');
        rois_sort=names(INDSORT);
       
        set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',rois_sort, 'XTick',1:numel(rois_sort));
          rotateXLabels(gca(), 90 );
           set(gca, 'YTickLabel',rois_sort, 'YTick',1:numel(rois_sort));
        hold off;                         
        
        %% %% Participation Coefficient        
        [Ppos,Pneg]=participation_coef_sign(fc,Ci0);
        [Ppos,I]=sort(Ppos,1,'descend'); 
	       ROIs_sorted=names(I);
           figure;
           bar(Ppos);
          title('Positive Clustering Coefficient');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
          
          [Pneg,I]=sort(Pneg,1,'descend'); 
	       ROIs_sorted=names(I);
           figure;
           bar(Pneg);
          title('Negative Clustering Coefficient');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
    
 %% Effect of thresholding & binarizing
 
%       for i=0.01:0.01:0.1
%         th_fc=threshold_proportional(fc,i);
%         bin_fc=weight_conversion(th_fc,'binarize');
%         d=sum(bin_fc);
%         histogram(d,10); title([i]);
%         
%         pause(0.5);
%       end  
%       
      th_fc=threshold_proportional(fc,0.25);
      bin_fc=weight_conversion(th_fc,'binarize');
figure;
imagesc(th_fc);
figure;
imagesc(bin_fc);
 %% Centrality Measures %%
      % Betweenness Centrality %
           b=betweenness_bin(bin_fc);
           norm_b=b./((N-1)*(N-2));
	       [B,I]=sort(norm_b,2,'descend'); 
	       ROIs_sorted=names(I);
           Bplot= bar(B);
          title('Betweenness Centrality');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
%           rotateXLabels(gca(), 90 );
         % setbarcolor(Bplot,B,y);           


      % Eigen Centrality %
          eig=eigenvector_centrality_und(th_fc);
          [E,I]=sort(eig,1,'descend'); 
	      ROIs_sorted=names(I);
          figure;
          Eplot= bar(E);
          title('Eigenvector Centrality');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
%           rotateXLabels( gca(), 90 );

%%%% See tutorial for visualization of network %%%%
  %% Length Matrix - inverse of adj weights
           L=weight_conversion(th_fc,'lengths');
        
  %% Distance Matrix - Shortest path for each node pair
           D= distance_bin(bin_fc);
        
  %% Efficiency 
          Eglob=efficiency_bin(bin_fc);
          Eloc= efficiency_bin(bin_fc,1);
          [SEloc,I]=sort(Eloc,1,'descend'); 
	       ROIs_sorted=names(I);
           Lplot= bar(SEloc);
          title('Local Efficiency');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
         
  %% Small World Index
         C=   clustering_coef_bu(bin_fc); %Clustering Coefficient
         norm_CC=C./((N-1)*(N-2));
         cpl= charpath(D,0,0); % Characteristic Path Length
         small_worldness= C./cpl;
         kden=density_und(bin_fc);
         %deg=Degrees_und   
 %% save Th_fc            
 %fid=fopen('TH_FC.edge','wt');
 writematrix(th_fc,"FC.txt");
 