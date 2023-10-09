function(add_doc_deployment target)
  set(options)
  set(one_value_args PAGES_REPO)
  set(multi_value_args)

  cmake_parse_arguments(deploy "${options}" "${one_value_args}"
    "${multi_value_args}" ${ARGN})

  if(NOT DEFINED deploy_PAGES_REPO)
    message(FATAL_ERROR "Need to specify PAGES_REPO")
  endif()

  find_package(Git REQUIRED)
  include(FleCSI/colors)

  _flecsi_define_doc_group_target()

  #--------------------------------------------------------------------------#
  # This target will work with multiple doxygen targets. However, because
  # sphinx is used for the main html page content, it will only work with
  # one sphinx target, i.e., the one named `sphinx`.
  #--------------------------------------------------------------------------#

  add_custom_target(${target}
    COMMAND
      echo "${FLECSI_Green}Updating pages branch${FLECSI_ColorReset}" &&
        ([ -e pages ] || (echo "cloning ${deploy_PAGES_REPO}" &&
          ${GIT_EXECUTABLE} clone -q --single-branch --branch pages
            ${deploy_PAGES_REPO} pages &&
          cd pages &&
          ${GIT_EXECUTABLE} rm -qr . && ${GIT_EXECUTABLE} reset -q &&
          ${GIT_EXECUTABLE} checkout .gitignore &&
          ${GIT_EXECUTABLE} checkout .gitlab-ci.yaml)) &&
        echo "${FLECSI_Green}Copying pages${FLECSI_ColorReset}" &&
        mkdir -p pages/public &&
        cp -rT doc pages/public &&
        (cd pages && ${GIT_EXECUTABLE} add -A .) &&
        echo "${FLECSI_Yellow}Updated pages are in"
          "${CMAKE_BINARY_DIR}/pages${FLECSI_ColorReset}"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  add_dependencies(${target} doc)
endfunction()
