function [sorted_string_out] = sort_complex_names(unsorted_string_in)
% Takes in a string including complexes in unsorted form (e.g. [B_C_A])
% Returns a string describing with complex-compounds in alphabetically sorted form (e.g. [A_B_C]) 

% Split the input string into an array of strings: [reactants; a1b1; complex; k1; product] 
unsorted_array=split(unsorted_string_in);

% The sorted array will differ from the unsorted array on rows 1,4,6 ([reactants, complex, product])
sorted_array = unsorted_array;
sorted_array(1)="";
sorted_array(4)="";
sorted_array(6)="";

%for reactants, complex, product:
for row=1:6
    if(row==1 || row==4 || row==6)
    unsorted_row_string = string(unsorted_array(row));
  
    %split the reactants/complex/product into all its compounds
    unsorted_row_string_array = split(unsorted_row_string,"+");

    %for every compound j in the reactants/complex/product
    for j = 1:size(unsorted_row_string_array,1)
        %define
        unsorted_row_string_array(j);
        %split and sort
        sorted_compound_array = extractBetween((unsorted_row_string_array(j)),"[","]");
        sorted_compound_array = split(sorted_compound_array,'_');
        sorted_compound_array = sort(sorted_compound_array);
        %reconstruct
        sorted_compound_string = "[";
        for jj = 1:size(sorted_compound_array, 1)
            sorted_compound_string = sorted_compound_string + sorted_compound_array(jj);
            if(jj<size(sorted_compound_array, 1))
               sorted_compound_string = sorted_compound_string + "_";
            end
        end
        sorted_compound_string = sorted_compound_string + "]";
        sorted_array(row)=sorted_array(row)+sorted_compound_string;
        if(j<size(unsorted_row_string_array,1))
            sorted_array(row)=sorted_array(row)+'+';
        end
    end 
    end
  
end
%join the sorted array to one output string
sorted_string_out=strjoin(sorted_array);
end

