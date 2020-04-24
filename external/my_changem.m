function mapout = my_changem(mapout, newcode, oldcode)
   assert(numel(newcode) == numel(oldcode), 'newcode and oldecode must have the same number of elements');
   [toreplace, bywhat] = ismember(mapout, oldcode);
   mapout(toreplace) = newcode(bywhat(toreplace));
end