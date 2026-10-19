emcc ecm.cpp mini-gmp.c -O3 \
  --bind \
  -o ecm_module.js \
  -s WASM=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s SINGLE_FILE=1
