cmake_minimum_required(VERSION 3.18)
project(cudaSAXS_GUI)

# Find Qt6 or Qt5
find_package(Qt6 COMPONENTS Core Widgets QUIET)
if (NOT Qt6_FOUND)
    find_package(Qt5 5.15 REQUIRED COMPONENTS Core Widgets)
    set(QT_VERSION_MAJOR 5)
else()
    set(QT_VERSION_MAJOR 6)
endif()

# Set up Qt
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Include directories
include_directories(include)

# Source files
set(SOURCES
    src/main.cpp
    src/MainWindow.cpp
    src/InputForm.cpp
    src/AdvancedOptionsDialog.cpp
    src/ProcessRunner.cpp
    src/OptionsData.cpp
)

# Header files
set(HEADERS
    include/MainWindow.h
    include/InputForm.h
    include/AdvancedOptionsDialog.h
    include/ProcessRunner.h
    include/OptionsData.h
)

# Create executable
add_executable(cudaSAXS_GUI ${SOURCES} ${HEADERS})

# Link Qt libraries
if(QT_VERSION_MAJOR EQUAL 6)
    target_link_libraries(cudaSAXS_GUI Qt6::Core Qt6::Widgets)
else()
    target_link_libraries(cudaSAXS_GUI Qt5::Core Qt5::Widgets)
endif()

# Installation
install(TARGETS cudaSAXS_GUI DESTINATION bin)

# On Windows, deploy Qt libraries
if(WIN32)
    if(QT_VERSION_MAJOR EQUAL 6)
        add_custom_command(TARGET cudaSAXS_GUI POST_BUILD
            COMMAND ${CMAKE_PREFIX_PATH}/bin/windeployqt.exe $<TARGET_FILE:cudaSAXS_GUI>
        )
    else()
        add_custom_command(TARGET cudaSAXS_GUI POST_BUILD
            COMMAND ${CMAKE_PREFIX_PATH}/bin/windeployqt.exe $<TARGET_FILE:cudaSAXS_GUI>
        )
    endif()
endif()