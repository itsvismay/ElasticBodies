FILE(REMOVE_RECURSE
  "libglfw3.pdb"
  "libglfw3.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/glfw.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
