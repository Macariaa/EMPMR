function crosssection( MShape,ind,Vmin,Vmax )
id= [1,2,3];
id(ind)= [];
figure 
[va,idd]= find(MShape(ind,:)<Vmax);
MData= MShape(:,idd);
[va,idd]= find(MData(ind,:)>Vmin);
MData= MData(:,idd);
plot(MData(id(2),:),MData(id(1),:),'.b')
end

