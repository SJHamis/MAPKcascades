function [conservation_laws_math_LHS_out, conservation_laws_math_RHS_out, conservation_laws_text_out, super_compound_list_index]...
    =  get_conservation_laws(list_of_super_compounds_in, list_of_sub_compounds_in)

   %Remove the extra slashes \[A\] becomes [A]
   list_of_sub_compounds_in = regexprep(list_of_sub_compounds_in, '\\[', '[');
   list_of_sub_compounds_in = regexprep(list_of_sub_compounds_in, '\\]', ']');
   
   %Supercompounds are in form "A"
   
   conservation_laws_text_out="";
   %Loop throug all super compounds
   for i = 1:size(list_of_super_compounds_in,1)
       temp_string="";
       temp_super_compound=list_of_super_compounds_in(i);
       %Loop throug all sub compounds
       for j = 1:size(list_of_sub_compounds_in,1)
           temp_sub_compound=list_of_sub_compounds_in(j);
           %if super appears in sub
           if(contains(temp_sub_compound, temp_super_compound ))
               if(strcmp(temp_sub_compound,"["+temp_super_compound+"]")==0)%if it is not itself
                    %temp_string=temp_string + "-" + temp_sub_compound;
                    temp_string=temp_string + "+" + temp_sub_compound;
               end
           end
       end
       
       %temp_string= "[" + temp_super_compound + "]=" + temp_super_compound+"_tot" + temp_string;
       temp_string= "0 = [" + temp_super_compound + "]" + temp_string + "-" + temp_super_compound+"_tot" ;
       
       conservation_laws_text_out=[conservation_laws_text_out;temp_string];
   end  %End of loop throug all super compounds
   
   %removes the top (empty) line
   conservation_laws_text_out(cellfun('isempty',conservation_laws_text_out))=[];
   
   conservation_laws_math_RHS_out =  conservation_laws_text_out;
   %%Mathify! for all of the compounds
   for i = 1:length(list_of_sub_compounds_in)
        temp_string=list_of_sub_compounds_in(i);
        i_string=string(i);
        conservation_laws_math_RHS_out = strrep(conservation_laws_math_RHS_out, temp_string,  "y("+i_string+")");
   end

  conservation_laws_math_RHS_out = strrep(conservation_laws_math_RHS_out, "0 = ","");
  conservation_laws_math_LHS_out = string(zeros(size(conservation_laws_math_RHS_out,1),1));
  
  super_compound_list_index=zeros(1,size(list_of_super_compounds_in,1));
  index=1;
  for j = 1:size(list_of_super_compounds_in,1)
        for i = 1:size(list_of_sub_compounds_in,1)      
           if( strcmp("["+list_of_super_compounds_in(j)+"]",list_of_sub_compounds_in(i)) )
               super_compound_list_index(index)=i;
               index=index+1;
           end
       end
  end 
  
  super_compound_list_index;
end

   
