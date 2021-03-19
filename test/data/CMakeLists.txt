# TODO - create data.csv if not present
configure_file(data.csv ${CMAKE_BINARY_DIR}/bin/data.csv COPYONLY)
configure_file(test.bgen ${CMAKE_BINARY_DIR}/bin/test.bgen COPYONLY)
configure_file(test.bgen.bgi ${CMAKE_BINARY_DIR}/bin/test.bgen.bgi COPYONLY)
