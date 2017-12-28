# This CMake script will executed to generate the project specific
# signature files in case no version control is present
# 
# @author: Norbert Manthey <nmanthey@conp-solutions.com>

set(TEMP_FILE ${DESTINATION}.tmp)

# Let the project-specific CMake script build the version file
execute_process(COMMAND ${CMAKE_COMMAND} -D PROJECT_NAME=${PROJECT_NAME}
                                         -D OUTPUT=${TEMP_FILE}
                                         -D GIT_SHA1=unavailable
                                         -D GIT_DATE=unavailable
                                         -D VERSION=${VERSION}
                                         -P ${SCRIPT})
# copy the generated version file only if the commit (Git SHA1) has changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
                                         ${TEMP_FILE}
                                         ${DESTINATION})
# Remove the temporary generated file
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${TEMP_FILE})
