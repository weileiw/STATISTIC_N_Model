function func = make_MLR_function(mdl,tbl)
func = 0;
VA = {'x1' 'x2' 'x3' 'x4' 'x5' 'x6'};
Var_names = mdl.CoefficientNames;
CoefficientValue = mdl.Coefficients.Estimate;

for jl = 1:length(Var_names)
    if strcmp(Var_names{jl},'(Intercept)') 
        func = func + CoefficientValue(jl);
       
    end
end 

for j1 = 1:length(Var_names)
    for j2 = 1:length(VA)
        if strcmp(Var_names{j1},VA{j2})
            func = func + CoefficientValue(j1)*tbl(:,j2);
            
        end
    end
end

for jl = 1:length(Var_names)
    if strfind(Var_names{jl},':')>0
        tmp = Var_names{jl};
        for jk = 1:length(VA)
            if  strcmp(VA{jk},tmp(1:2));
                data_row1 = jk;
            end 
            if  strcmp(VA{jk},tmp(4:5));
                data_row2 = jk;
            end     
        end
        func = func + CoefficientValue(jl)* tbl(:,data_row1).*tbl(:,data_row2);
        
    end
end 

for jl = 1:length(Var_names)
    if strfind(Var_names{jl},'^2')
        tmp = Var_names{jl};
        for jm = 1:length(VA)
            if strcmp(VA{jm},tmp(1:2));
                func = func + CoefficientValue(jl)*tbl(:,jm).^2;
      
            end
        end
    end 
end

