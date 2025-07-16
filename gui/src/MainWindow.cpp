#include "MainWindow.h"
#include "InputForm.h"
#include "ProcessRunner.h"
#include "AdvancedOptionsDialog.h"
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QToolBar>
#include <QStatusBar>
#include <QDockWidget>
#include <QTextEdit>
#include <QProgressBar>
#include <QLabel>
#include <QMessageBox>
#include <QFileDialog>
#include <QSettings>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QFile>
#include <QDateTime>
#include <QDir>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
      m_processRunner(new ProcessRunner(this)),
      m_advancedDialog(new AdvancedOptionsDialog(this)),
      m_elapsedTimer(new QTimer(this)),
      m_startTime(0)
{

    setWindowTitle("cudaSAXS GUI");
    resize(800, 600);

    setupUI();
    createMenus();
    createStatusBar();
    connectSignals();

    // Load settings
    QSettings settings("cudaSAXS", "GUI");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());

    // Set executable path (adjust as needed)
    QString execPath = settings.value("executablePath", "./cudaSAXS").toString();
    m_processRunner->setExecutablePath(execPath);
}

MainWindow::~MainWindow()
{
    // Save settings
    QSettings settings("cudaSAXS", "GUI");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
}

void MainWindow::setupUI()
{
    // Central widget
    m_inputForm = new InputForm(this);
    setCentralWidget(m_inputForm);

    // Output dock
    auto outputDock = new QDockWidget("Output", this);
    m_outputText = new QTextEdit;
    m_outputText->setReadOnly(true);
    m_outputText->setFont(QFont("Consolas", 9));
    outputDock->setWidget(m_outputText);
    addDockWidget(Qt::BottomDockWidgetArea, outputDock);
}

void MainWindow::createMenus()
{
    // File menu
    auto fileMenu = menuBar()->addMenu("&File");

    auto newAction = fileMenu->addAction("&New Project");
    newAction->setShortcut(QKeySequence::New);
    connect(newAction, &QAction::triggered, this, &MainWindow::onNewProject);

    auto openAction = fileMenu->addAction("&Open Project...");
    openAction->setShortcut(QKeySequence::Open);
    connect(openAction, &QAction::triggered, this, &MainWindow::onOpenProject);

    fileMenu->addSeparator();

    m_saveAction = fileMenu->addAction("&Save Project");
    m_saveAction->setShortcut(QKeySequence::Save);
    connect(m_saveAction, &QAction::triggered, this, &MainWindow::onSaveProject);

    auto saveAsAction = fileMenu->addAction("Save Project &As...");
    saveAsAction->setShortcut(QKeySequence::SaveAs);
    connect(saveAsAction, &QAction::triggered, this, &MainWindow::onSaveProjectAs);

    fileMenu->addSeparator();

    auto exitAction = fileMenu->addAction("E&xit");
    exitAction->setShortcut(QKeySequence::Quit);
    connect(exitAction, &QAction::triggered, this, &QWidget::close);

    // Run menu
    auto runMenu = menuBar()->addMenu("&Run");

    m_runAction = runMenu->addAction("&Run SAXS");
    m_runAction->setShortcut(Qt::Key_F5);
    connect(m_runAction, &QAction::triggered, this, &MainWindow::onRunRequested);

    m_stopAction = runMenu->addAction("&Stop");
    m_stopAction->setShortcut(Qt::Key_F6);
    m_stopAction->setEnabled(false);
    connect(m_stopAction, &QAction::triggered, m_processRunner, &ProcessRunner::cancel);

    // Tools menu
    auto toolsMenu = menuBar()->addMenu("&Tools");

    auto advancedAction = toolsMenu->addAction("&Advanced Options...");
    connect(advancedAction, &QAction::triggered, this, &MainWindow::onAdvancedOptionsRequested);

    toolsMenu->addSeparator();

    auto prefsAction = toolsMenu->addAction("&Preferences...");
    connect(prefsAction, &QAction::triggered, this, &MainWindow::onPreferences);

    // Help menu
    auto helpMenu = menuBar()->addMenu("&Help");

    auto aboutAction = helpMenu->addAction("&About...");
    connect(aboutAction, &QAction::triggered, this, &MainWindow::onAbout);

    // Toolbar
    auto toolbar = addToolBar("Main");
    toolbar->addAction(m_runAction);
    toolbar->addAction(m_stopAction);
    toolbar->addSeparator();
    toolbar->addAction(newAction);
    toolbar->addAction(openAction);
    toolbar->addAction(m_saveAction);
}

void MainWindow::createStatusBar()
{
    m_statusLabel = new QLabel("Ready");
    statusBar()->addWidget(m_statusLabel);

    m_progressBar = new QProgressBar;
    m_progressBar->setVisible(false);
    m_progressBar->setTextVisible(false);
    m_progressBar->setRange(0, 0); // Indeterminate
    statusBar()->addPermanentWidget(m_progressBar);

    m_elapsedLabel = new QLabel;
    statusBar()->addPermanentWidget(m_elapsedLabel);
}

void MainWindow::connectSignals()
{
    // Input form signals
    connect(m_inputForm, &InputForm::runRequested,
            this, &MainWindow::onRunRequested);
    connect(m_inputForm, &InputForm::advancedOptionsRequested,
            this, &MainWindow::onAdvancedOptionsRequested);
    connect(m_inputForm, &InputForm::formValidityChanged,
            this, &MainWindow::onFormValidityChanged);

    // Process runner signals
    connect(m_processRunner, &ProcessRunner::started,
            this, &MainWindow::onProcessStarted);
    connect(m_processRunner, &ProcessRunner::finished,
            this, &MainWindow::onProcessFinished);
    connect(m_processRunner, &ProcessRunner::error,
            this, &MainWindow::onProcessError);
    connect(m_processRunner, &ProcessRunner::outputReceived,
            this, &MainWindow::onProcessOutput);
    connect(m_processRunner, &ProcessRunner::progressReceived,
            this, &MainWindow::onProcessProgress);

    // Timer
    connect(m_elapsedTimer, &QTimer::timeout,
            this, &MainWindow::updateElapsedTime);
}

void MainWindow::onRunRequested()
{
    if (m_processRunner->isRunning())
    {
        QMessageBox::warning(this, "Process Running",
                             "A calculation is already running. Please wait for it to finish or stop it.");
        return;
    }

    OptionsData options = m_inputForm->getOptionsData();
    if (!options.isValid())
    {
        QMessageBox::critical(this, "Invalid Options", options.validationError());
        return;
    }

    // Clear output
    m_outputText->clear();

    // Get command line arguments
    QStringList args = options.toCommandLineArgs();

    // Log the command
    m_outputText->append("Command: cudaSAXS " + args.join(" ") + "\n");
    m_outputText->append(QString(80, '-') + "\n");

    // Set working directory to the directory of the trajectory file
    QFileInfo xtcInfo(options.xtcFile);
    m_processRunner->setWorkingDirectory(xtcInfo.absolutePath());

    // Run the process
    m_processRunner->run(args);
}

void MainWindow::onAdvancedOptionsRequested()
{
    OptionsData currentOptions = m_inputForm->getOptionsData();
    m_advancedDialog->setOptionsData(currentOptions);

    if (m_advancedDialog->exec() == QDialog::Accepted)
    {
        OptionsData advancedOptions = m_advancedDialog->getOptionsData();
        // Merge with current options
        currentOptions.bsplineOrder = advancedOptions.bsplineOrder;
        currentOptions.scaleFactor = advancedOptions.scaleFactor;
        currentOptions.scaledGrid = advancedOptions.scaledGrid;
        currentOptions.waterModel = advancedOptions.waterModel;
        currentOptions.sodiumAtoms = advancedOptions.sodiumAtoms;
        currentOptions.chlorineAtoms = advancedOptions.chlorineAtoms;

        m_inputForm->setOptionsData(currentOptions);
    }
}

void MainWindow::onProcessStarted()
{
    m_runAction->setEnabled(false);
    m_stopAction->setEnabled(true);
    m_progressBar->setVisible(true);
    m_statusLabel->setText("Running...");

    // Start elapsed timer
    m_startTime = QDateTime::currentMSecsSinceEpoch();
    m_elapsedTimer->start(1000); // Update every second
    updateElapsedTime();
}

void MainWindow::onProcessFinished(int exitCode)
{
    m_runAction->setEnabled(true);
    m_stopAction->setEnabled(false);
    m_progressBar->setVisible(false);
    m_elapsedTimer->stop();

    if (exitCode == 0)
    {
        m_statusLabel->setText("Finished successfully");
        m_outputText->append("\n" + QString(80, '-'));
        m_outputText->append("Process finished successfully");

        QMessageBox::information(this, "Success",
                                 "SAXS calculation completed successfully!");
    }
    else
    {
        m_statusLabel->setText("Finished with errors");
        m_outputText->append("\n" + QString(80, '-'));
        m_outputText->append(QString("Process finished with exit code: %1").arg(exitCode));

        QMessageBox::warning(this, "Process Failed",
                             QString("Process exited with code %1. Check the output for details.").arg(exitCode));
    }
}

void MainWindow::onProcessError(const QString &error)
{
    m_outputText->append("\nERROR: " + error);
    m_statusLabel->setText("Error occurred");
}

void MainWindow::onProcessOutput(const QString &output)
{
    // Add output, but don't scroll if user is reading
    auto cursor = m_outputText->textCursor();
    cursor.movePosition(QTextCursor::End);
    cursor.insertText(output);
}

void MainWindow::onProcessProgress(const QString &progress)
{
    m_statusLabel->setText(progress);
}

void MainWindow::onFormValidityChanged(bool isValid)
{
    m_runAction->setEnabled(isValid && !m_processRunner->isRunning());
}

void MainWindow::updateElapsedTime()
{
    if (m_startTime > 0)
    {
        qint64 elapsed = QDateTime::currentMSecsSinceEpoch() - m_startTime;
        int seconds = elapsed / 1000;
        int minutes = seconds / 60;
        seconds %= 60;
        m_elapsedLabel->setText(QString("Elapsed: %1:%2")
                                    .arg(minutes, 2, 10, QChar('0'))
                                    .arg(seconds, 2, 10, QChar('0')));
    }
    else
    {
        m_elapsedLabel->clear();
    }
}

void MainWindow::onNewProject()
{
    if (m_processRunner->isRunning())
    {
        QMessageBox::warning(this, "Process Running",
                             "Cannot create new project while a calculation is running.");
        return;
    }

    m_projectFile.clear();
    m_inputForm->resetToDefaults();
    m_outputText->clear();
    setWindowTitle("cudaSAXS GUI");
}

void MainWindow::onOpenProject()
{
    if (m_processRunner->isRunning())
    {
        QMessageBox::warning(this, "Process Running",
                             "Cannot open project while a calculation is running.");
        return;
    }

    QString filename = QFileDialog::getOpenFileName(this,
                                                    "Open Project", QString(), "cudaSAXS Project (*.csp);;All Files (*)");

    if (!filename.isEmpty())
    {
        loadProject(filename);
    }
}

void MainWindow::onSaveProject()
{
    if (m_projectFile.isEmpty())
    {
        onSaveProjectAs();
    }
    else
    {
        saveProject(m_projectFile);
    }
}

void MainWindow::onSaveProjectAs()
{
    QString filename = QFileDialog::getSaveFileName(this,
                                                    "Save Project", QString(), "cudaSAXS Project (*.csp);;All Files (*)");

    if (!filename.isEmpty())
    {
        if (!filename.endsWith(".csp"))
        {
            filename += ".csp";
        }
        saveProject(filename);
    }
}

void MainWindow::loadProject(const QString &filename)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly))
    {
        QMessageBox::critical(this, "Error",
                              "Could not open project file: " + file.errorString());
        return;
    }

    QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
    QJsonObject root = doc.object();

    // Load options
    QJsonObject options = root["options"].toObject();
    OptionsData data;

    data.tprFile = options["tprFile"].toString();
    data.xtcFile = options["xtcFile"].toString();

    QJsonArray gridArray = options["gridSize"].toArray();
    for (const auto &val : gridArray)
    {
        data.gridSize.append(val.toInt());
    }

    data.startFrame = options["startFrame"].toInt();
    data.endFrame = options["endFrame"].toInt();
    data.frameInterval = options["frameInterval"].toInt(1);
    data.outputFile = options["outputFile"].toString("saxs.dat");
    data.bsplineOrder = options["bsplineOrder"].toInt(4);

    QJsonArray scaledArray = options["scaledGrid"].toArray();
    for (const auto &val : scaledArray)
    {
        data.scaledGrid.append(val.toInt());
    }

    data.scaleFactor = options["scaleFactor"].toDouble(2.5);
    data.binSize = options["binSize"].toDouble(0.05);
    data.qCutoff = options["qCutoff"].toDouble(4.0);
    data.waterModel = options["waterModel"].toString();
    data.sodiumAtoms = options["sodiumAtoms"].toInt();
    data.chlorineAtoms = options["chlorineAtoms"].toInt();
    data.simulationType = options["simulationType"].toString("npt");

    m_inputForm->setOptionsData(data);

    // Load output if present
    if (root.contains("output"))
    {
        m_outputText->setPlainText(root["output"].toString());
    }

    m_projectFile = filename;
    setWindowTitle(QString("cudaSAXS GUI - %1").arg(QFileInfo(filename).fileName()));
}

void MainWindow::saveProject(const QString &filename)
{
    QJsonObject root;

    // Save options
    OptionsData data = m_inputForm->getOptionsData();
    QJsonObject options;

    options["tprFile"] = data.tprFile;
    options["xtcFile"] = data.xtcFile;

    QJsonArray gridArray;
    for (int val : data.gridSize)
    {
        gridArray.append(val);
    }
    options["gridSize"] = gridArray;

    options["startFrame"] = data.startFrame;
    options["endFrame"] = data.endFrame;
    options["frameInterval"] = data.frameInterval;
    options["outputFile"] = data.outputFile;
    options["bsplineOrder"] = data.bsplineOrder;

    QJsonArray scaledArray;
    for (int val : data.scaledGrid)
    {
        scaledArray.append(val);
    }
    options["scaledGrid"] = scaledArray;

    options["scaleFactor"] = data.scaleFactor;
    options["binSize"] = data.binSize;
    options["qCutoff"] = data.qCutoff;
    options["waterModel"] = data.waterModel;
    options["sodiumAtoms"] = data.sodiumAtoms;
    options["chlorineAtoms"] = data.chlorineAtoms;
    options["simulationType"] = data.simulationType;

    root["options"] = options;

    // Save output
    root["output"] = m_outputText->toPlainText();

    // Write to file
    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly))
    {
        QMessageBox::critical(this, "Error",
                              "Could not save project file: " + file.errorString());
        return;
    }

    QJsonDocument doc(root);
    file.write(doc.toJson());

    m_projectFile = filename;
    setWindowTitle(QString("cudaSAXS GUI - %1").arg(QFileInfo(filename).fileName()));
    m_statusLabel->setText("Project saved");
}

void MainWindow::onAbout()
{
    QMessageBox::about(this, "About cudaSAXS GUI",
                       "<h3>cudaSAXS GUI</h3>"
                       "<p>A graphical interface for cudaSAXS - GPU-accelerated SAXS calculations</p>"
                       "<p>This interface provides an easy way to configure and run SAXS calculations "
                       "on molecular dynamics trajectories using CUDA acceleration.</p>"
                       "<p>Version 1.0</p>");
}

void MainWindow::onPreferences()
{
    QSettings settings("cudaSAXS", "GUI");

    QString currentPath = settings.value("executablePath", "./cudaSAXS").toString();
    QString newPath = QFileDialog::getOpenFileName(this,
                                                   "Select cudaSAXS Executable", currentPath, "Executable Files (*)");

    if (!newPath.isEmpty())
    {
        settings.setValue("executablePath", newPath);
        m_processRunner->setExecutablePath(newPath);
        QMessageBox::information(this, "Preferences",
                                 "Executable path updated. The new path will be used for future calculations.");
    }
}