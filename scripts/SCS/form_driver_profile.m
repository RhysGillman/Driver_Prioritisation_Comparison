function [ driver_module ] = form_driver_profile( driver_control_module,predict_driver )
%只保留driver_control_module中和predict_driver有交集的行
%Only keep the line in both driver_control_module and predict_driver 
driver_module=[];
[row_driver_control_module,~]=size(driver_control_module);
k=1;
for i=1:row_driver_control_module
    
   if ismember(driver_control_module{i,1},predict_driver)~=0
       driver_module{k,1}=driver_control_module{i,1};
       driver_module{k,2}=driver_control_module{i,2};
       driver_module{k,3}=driver_control_module{i,3};
       k=k+1;
   end
    
end


end

