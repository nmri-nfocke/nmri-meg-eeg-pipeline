function txt=legalize_label(txt)
 txt=strtok(txt,'(');
 txt=strtrim(txt);
 txt=strrep(txt,'>','more');
 txt=strrep(txt,'<','less'); 
 txt=strrep(txt,' - ','-');
 txt=strrep(txt,'- ','-');
 txt=strrep(txt,' -','-');
 txt=strrep(txt,' ','_');
 txt=strrep(txt,'[','_');
 txt=strrep(txt,']','');
 
 txt=strrep(txt,'__','_');
 if strcmp(txt(end),'_')
  txt=txt(1:end-1);
 end
end
