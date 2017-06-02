FILE(REMOVE_RECURSE
  "libnanogui.pdb"
  "libnanogui.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/nanogui.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
