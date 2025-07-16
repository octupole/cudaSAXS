#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include "OptionsData.h"

QT_BEGIN_NAMESPACE
class QTextEdit;
class QProgressBar;
class QLabel;
class QAction;
QT_END_NAMESPACE

class InputForm;
class ProcessRunner;
class AdvancedOptionsDialog;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void onRunRequested();
    void onAdvancedOptionsRequested();
    void onProcessStarted();
    void onProcessFinished(int exitCode);
    void onProcessError(const QString &error);
    void onProcessOutput(const QString &output);
    void onProcessProgress(const QString &progress);
    void onFormValidityChanged(bool isValid);

    // Menu actions
    void onNewProject();
    void onOpenProject();
    void onSaveProject();
    void onSaveProjectAs();
    void onAbout();
    void onPreferences();

    // Timer for status
    void updateElapsedTime();

private:
    void setupUI();
    void createMenus();
    void createStatusBar();
    void connectSignals();

    void loadProject(const QString &filename);
    void saveProject(const QString &filename);

    // Central widget
    InputForm *m_inputForm;

    // Output dock
    QTextEdit *m_outputText;

    // Process runner
    ProcessRunner *m_processRunner;

    // Advanced options dialog
    AdvancedOptionsDialog *m_advancedDialog;

    // Status bar widgets
    QLabel *m_statusLabel;
    QProgressBar *m_progressBar;
    QLabel *m_elapsedLabel;

    // Timer for elapsed time
    QTimer *m_elapsedTimer;
    qint64 m_startTime;

    // Current project file
    QString m_projectFile;

    // Menu actions
    QAction *m_runAction;
    QAction *m_stopAction;
    QAction *m_saveAction;
};

#endif // MAINWINDOW_H